"""Approximation error of a mixed-logit model to rho^pi, by EM (fish, no fixed effects).

Usage:
    python approx_rum_EM_vertices_fish.py <index> <order> [data.csv]

    index   0-23, position in permutations((0, 1, 2, 3))
    order   1 for linear in characteristics, 2 for quadratic

Writes OUT_DIR/Fish_em_<order>_<index>.txt holding pi(0..3) followed by the
approximation error, minimised over EM_STARTS random initial points.

Environment: OUT_DIR, FISH_DATA, SEED_BASE, EM_STARTS, EM_MIXTURE, EM_MAX_ITER,
EM_TOL, EM_STANDARDIZE.
"""
import numpy as np
import pandas as pd
from itertools import combinations, permutations
from scipy.optimize import minimize, root
from tqdm import tqdm
from joblib import Parallel, delayed
from collections import defaultdict
from collections.abc import Iterable
import sys
import time
import warnings
import random
warnings.simplefilter('ignore')

index = int(sys.argv[1])
order = int(sys.argv[2])

import os
def _load_features(default='../../data/fish/fish_J4_K2.csv', argv_pos=3):
    path = sys.argv[argv_pos] if len(sys.argv) > argv_pos else os.environ.get('FISH_DATA', default)
    df = pd.read_csv(path)
    print('characteristics from', path, flush=True)
    print(df.to_string(index=False), flush=True)
    return df.drop(columns=['alt_id']).to_numpy(dtype=float)

def _out_dir():
    d = os.environ.get('OUT_DIR', 'output/raw')
    os.makedirs(d, exist_ok=True)
    return d

def _n_starts(default=5):
    """Random initialisations of EM per ranking; the reported error is the minimum."""
    return int(os.environ.get('EM_STARTS', default))

def _tol(default=1e-6):
    """Termination threshold on the absolute change in error between EM iterations."""
    return float(os.environ.get('EM_TOL', default))

def _standardise():
    """Whether to centre and scale the feature columns after any order-2 expansion."""
    return os.environ.get('EM_STANDARDIZE', '0') not in ('0', '', 'false', 'False')

def _max_iter(default=50):
    """Iteration cap per EM fit."""
    return int(os.environ.get('EM_MAX_ITER', default))

def _lambda_floor(default=0.0):
    """Lower bound on each mixture weight, renormalised. 0 disables it."""
    return float(os.environ.get('EM_LAMBDA_FLOOR', default))

def _mixture(default=10):
    """Number of mixture components M."""
    return int(os.environ.get('EM_MIXTURE', default))

_SEED = int(os.environ.get('SEED_BASE', 20240808)) + 1000 * order + index
np.random.seed(_SEED)
random.seed(_SEED)
print('seed', _SEED, flush=True)

#print(index)

#fixed_effects_dict = [[x,y,z] for x in range(-10,11) for y in range(-10,11) for z in range(-10,11) ]
#fixed_effects= np.array(fixed_effects_dict[index])

#fixed_effects = np.append(fixed_effects,0)
#fixed_effects = np.array(fixed_effects)
class ApproxRUM():
    """
    Attributes
    ----------
        K : int
            dimension of features
        X : tuple
            like (0, 1, ..., K-1)
        cal_D : list of tuples
            list of choice sets
    """
    def __init__(self, features):
        """
        Parameters
        ----------
            features : array
                features of each alternative
        """
        self.features = features
        self.K = self.features.shape[1]
        self.X = tuple(range(len(self.features)))
        self.cal_D = list(combinations(self.X, 2)) + [self.X] + list(combinations(self.X, 3))
        self.num_choice_sets = len(self.cal_D)

    def logit(self, D, beta, fixed_effects=0):
        """
        choice probabilities based on a logit model

        Parameters
        ----------
            D : tuple
                a choice set that contains X_ALT
            beta : K-dim array
                coefficients of logit model
            fixed_effects : len(D)-dim array, default 0
                fixed effects

        Returns
        -------
            probs : array
                an np.array of probabilities with length len(D)
        """
        num_alts = len(D)
        K = len(beta)
        D_feat = self.features[tuple([D])]
        assert D_feat.shape == (num_alts, K)
        lin_part_vec = D_feat @ beta + fixed_effects
        lin_part_vec_normalized = lin_part_vec - max(lin_part_vec)
        pre_probs = np.exp(lin_part_vec_normalized)
        probs = pre_probs / pre_probs.sum()
        return probs

    def mixed_logit(self, D, beta_mat, lambda_vec, fixed_effects=0):
        """
        choice probabilities based on a mixed logit model

        Parameters
        ----------
            D : tuple
                choice set
            beta_mat : (M, K)-dim array
                matrix of coeeficients
            lambda_vec : M-dim array
                M-dimensional vector of mixture weights
            fixed_effects : len(D)-dim array, default 0
                fixed effects

        Returns
        -------
            probs : array
                an np.array of probabilities with length len(D)
        """
        probs = 0
        for m, lam in enumerate(lambda_vec):
            beta = beta_mat[m, :]
            probs += lam * self.logit(D, beta, fixed_effects)
        return probs

    def get_gamma(self, D, lambda_old, beta_old_mat):
        """
        Parameters
        ----------
            lambda_old : M-dim array
                previous mixture weights
            beta_old_mat : (M, K)-dim array
                previous coefficients

        Returns
        -------
            gamma_mat : (M, len(D))-dim array
                an np.array of gammas with shape M \times len(D)
        """
        prob_vec = np.array([self.logit(D, beta) for beta in beta_old_mat])  # M \times len(D)
        weighted_prob_vec = prob_vec * lambda_old.reshape(-1, 1)  # M \times len(D)
        # A component whose weight underflows to zero can never recover: its next
        # lambda is built from a gamma row that is now identically zero. Once every
        # component in a column has collapsed this division is 0/0, and the NaN
        # propagates through lambda and every later gamma, so the whole fit returns
        # NaN. Guarding it leaves that column at zero -- finite, so the start simply
        # loses the minimum instead of poisoning it.
        denom = weighted_prob_vec.sum(axis=0)
        gamma_mat = weighted_prob_vec / np.where(denom > 0, denom, 1.0)
        return gamma_mat

    def get_lambda(self, gamma_new):
        """
        Parameters
        ----------
            gamma_new : list
                list of np.arrays in order of cal_D

        Returns
        -------
            lambda_vec : M-dim array
                an np.array of lambdas with length M
        """
        pre_lambda = sum([(gamma_new[i] * self.rho_set[D].reshape(1, -1)).sum(axis=1) for i, D in enumerate(self.cal_D)])
        lambda_vec = pre_lambda / self.num_choice_sets
        # Optional: keep every component alive. lambda sums to 1, so adding a floor
        # and renormalising stops a component dying permanently, which is what makes
        # M=18 collapse on targets the model can fit exactly. Off by default, since
        # switching it on perturbs every fit.
        floor = _lambda_floor()
        if floor > 0:
            lambda_vec = (lambda_vec + floor) / (1.0 + len(lambda_vec) * floor)
        return lambda_vec

    def _get_beta_mat(self, gamma_new, beta_initial_mat):
        """
        Find the optimal beta_mat naively; slow

        Parameters
        ----------
            gamma_new : list
                list of np.arrays in order of cal_D
            beta_initial_mat : (M, K)-dim array
                initial valus of beta

        Returns
        -------
            beta_new_mat : (M, K)-dim array
                an np.array with shape M \times K

        Notes
        -----
            This is no longer used. See get_beta_mat, which returns same outcomes but is much faster than this.
        """
        M = len(gamma_new[0])
        #fix m
        def ll(beta_m, m):
            ll_m = 0
            for i, D in enumerate(self.cal_D):
                gamma_new_m = gamma_new[i][m, :]  # len(D)
                lin_part_vec = self.features[tuple([D])] @ beta_m  # len(D)
                denom_logit = np.sum(np.exp(lin_part_vec))  # scalor
                log_denom_logit = np.log(denom_logit)  # scalor
                ll_m += sum(self.rho_set[D] * gamma_new_m * (lin_part_vec - log_denom_logit))
            return ll_m

        beta_new_mat = []
        ll_list = []
        for m in range(M):
            target_func = lambda beta_m: -ll(beta_m, m)
            # optimize target_func wrt beta_m
            method = "BFGS"
            res = minimize(target_func, beta_initial_mat[m, :], bounds=[(-20, 20) for _ in range(self.K)], method=method)
            # let beta_m_new the solution
            if res.success:
                beta_m_new = res.x
            else:
                beta_m_new = beta_initial_mat[m, :]
            beta_new_mat.append(beta_m_new)
            ll_list.append(res.fun)

        beta_new_mat = np.array(beta_new_mat)  # M \times K
        return beta_new_mat

    def get_beta_mat(self, gamma_new, beta_initial_mat):
        """
        Find the optimal beta_mat based on the FOC; fast

        Parameters
        ----------
            gamma_new : list
                list of np.arrays in order of cal_D
            beta_initial_mat : (M, K)-dim array
                initial valus of beta

        Returns
        -------
            beta_new_mat : (M, K)-dim array
                an np.array with shape M \times K
        """
        def eq(beta_m, m):
            total = 0
            for i, D in enumerate(self.cal_D):
                self.rho_set[D]
                gamma_new_m = gamma_new[i][m, :]
                nomin = (np.exp(self.features[tuple([D])] @ beta_m).reshape(-1, 1) * self.features[tuple([D])]).sum(axis=0)  # len(D) \times K
                denom = np.exp(self.features[tuple([D])] @ beta_m).sum()
                total += sum((self.rho_set[D] * gamma_new_m).reshape(-1, 1) * (self.features[tuple([D])] - nomin / denom))
            return total

        M = len(gamma_new[0])
        beta_new_mat = []
        for m in range(M):
            target_func = lambda beta_m: eq(beta_m, m)
            res = root(target_func, beta_initial_mat[m, :])
            if res.success:
                beta_m_new = res.x
            else:
                beta_m_new = beta_initial_mat[m, :]
            beta_new_mat.append(beta_m_new)

        beta_new_mat = np.array(beta_new_mat)  # M \times K
        return beta_new_mat

    def EM_algorithm(self, target_pref, M, max_iter_num,fixed_effects):
        """
        Parameters
        ----------
            target_pref : tuple
                a tuple with length n like (0, 1, 3, 2), which means 0 > 1 > 3 > 2
            M : int
                the number of mixtures
            max_iter_num : int
                number of iterarions of the EM algorithm

        Notes
        -----
            Fixed effects are not considered in this function.
        """
        # set the target preference
        self.set_rho_from_pref(target_pref)

        # initialize parameters
        lambda_old = np.zeros(M) + 1 / M
        beta_old_mat = np.random.normal(scale=0.5, size=(M, self.K))
        # place to store parameters
        self.gamma_list = [None]
        self.beta_mat_list = [beta_old_mat]
        self.lambda_list = [lambda_old]
        self.rho_list = [None]
        self.l2_error_list = [np.inf]
        
        for i in range(max_iter_num):
            # E step
            gamma_new = [self.get_gamma(D, lambda_old, beta_old_mat) for D in self.cal_D]  # len(cal_D) \times (M \times len(D))

            # M step
            lambda_new = self.get_lambda(gamma_new)
            beta_initial_mat = beta_old_mat
            beta_new_mat = self.get_beta_mat(gamma_new, beta_initial_mat)

            # progress just for log
            rho_est = self.get_rho_from_mixed_logit(lambda_new, beta_new_mat,fixed_effects)
            error = self.l2_distance(self.rho, rho_est)

            # store parameter records
            self.gamma_list.append(gamma_new)
            self.beta_mat_list.append(beta_new_mat)
            self.lambda_list.append(lambda_new)
            self.rho_list.append(rho_est)
            self.l2_error_list.append(error)

            # termination judgement
            if abs(error - self.l2_error_list[-2]) < _tol():
                # print(i)
                print('tolerance reached')
                break

            # update parameters
            lambda_old = lambda_new
            beta_old_mat = beta_new_mat

        # range(max_iter_num) leaves i at max_iter_num-1, so a capped-out fit is otherwise
        # indistinguishable from a converged one in the logs.
        if i == max_iter_num - 1:
            self.hit_cap = True
            print('Warning: EM did not converge (hit %d-iteration cap).' % max_iter_num,
                  flush=True)
        else:
            self.hit_cap = False

    def set_rho_from_pref(self, target_pref):
        """
        Parameters
        ----------
            target_pref : tuple
                a tuple with length n like (0, 1, 3, 2), which means 0 > 1 > 3 > 2
        """
        self.rho = dict()
        rank_array = np.argsort(target_pref)
        for D in self.cal_D:
            rho = {}
            min_rank = rank_array[tuple([D])].min()
            for x in D:
                self.rho[(x, D)] = 1 if rank_array[x] == min_rank else 0
        self.rho_set = dict([(D, np.array([self.rho[(x, D)] for x in D])) for D in self.cal_D])

    def l2_distance(self, rho_1, rho_2):
        """
        Parameters
        ----------
            rho_1 : dict
                a stochastic choice function
            rho_2 : dict
                a stochastic choice function

        Returns
        -------
            error : float
                the l2 distance between rho_1 and rho_2
        """
        count = 0
        squered_error_total = 0
        for D in self.cal_D:
           count += 1
           for x in D:
#                count += 1
                squered_error_total += (rho_1[(x, D)] - rho_2[(x, D)]) ** 2
        error = (squered_error_total ) ** 0.5/count
        return error

    def get_rho_from_logit(self, beta, fixed_effects):
        """
        Parameters
        ----------
            beta : K-dim array
                coefficients
            fixed_effects : len(X)-dim array, default 0
                fixed effects

        Returns
        -------
            rho_est : dict
                a stochastic choice function determined by a logit model with beta
        """
        #print(fixed_effects)
        if not isinstance(fixed_effects, Iterable):
            fixed_effects = np.zeros(len(self.X))
        rho_est = dict()
        for D in self.cal_D:
            probs = self.logit(D, beta, fixed_effects[[D]])
            for i, x in enumerate(D):
                rho_est[(x, D)] = probs[i]
        return rho_est

    def get_rho_from_mixed_logit(self, lambda_vec, beta_mat, fixed_effects):
        """
        Parameters
        ----------
            lambda_vec : M-dim array
                mixture weights
            beta_mat : (M, K)-dim array
                matrix of coefficients
            fixed_effects : len(X)-dim array, default 0
                fixed effects

        Returns
        -------
            rho_est : dict
                a stochastic choice function determined by a logit model with beta
        """
        #print(fixed_effects)
        if not isinstance(fixed_effects, Iterable):
            fixed_effects = np.zeros(len(self.X))
        rho_est = dict()
        for D in self.cal_D:
            probs = self.mixed_logit(D, beta_mat, lambda_vec, fixed_effects[[D]])
            for i, x in enumerate(D):
                rho_est[(x, D)] = probs[i]
        return rho_est

    def set_rho_from_pref_mixture(self, pref_1, pref_2, weight):
        """
        Parameters
        ----------
            pref_1 : tuple
                a tuple with length n like (0, 1, 3, 2), which means 0 > 1 > 3 > 2
            pref_2 : tuple
                a tuple with length n like (0, 1, 3, 2), which means 0 > 1 > 3 > 2
            weight : float
                weight for pref_1
        """
        self.rho = defaultdict(int)
        rank_array_1 = np.argsort(pref_1)
        rank_array_2 = np.argsort(pref_2)
        for D in self.cal_D:
            min_rank_1 = rank_array_1[tuple([D])].min()
            min_rank_2 = rank_array_2[tuple([D])].min()
            for x in D:
                if rank_array_1[x] == min_rank_1:
                    self.rho[(x, D)] += weight
                if rank_array_2[x] == min_rank_2:
                    self.rho[(x, D)] += 1 - weight
        self.rho_set = dict([(D, np.array([self.rho[(x, D)] for x in D])) for D in self.cal_D])

   
    def update_beta(self, target_rho_set, fixed_effects):
        """
        Parameters
        ----------
            target_rho_set : dict
                a dictionary of which keys are choice sets
            fixed_effects : len(X)-dim array, default 0
                fixed effects

        Returns
        -------
            beta_new : K-dim array

        """
        if not isinstance(fixed_effects, Iterable):
            fixed_effects = np.zeros(len(self.X))

        def eq(beta):
            return sum([np.square(target_rho_set[D] - self.logit(D, beta, fixed_effects[[D]])).sum() for D in self.cal_D])

        beta_init = np.random.normal(scale=0.1, size=self.K)
        method = "BFGS"
        res = minimize(eq, beta_init, bounds=[(-20, 20) for _ in range(self.K)], method=method)
        beta_new = res.x
        return beta_new


def process_EM(pref, features,fixed_effects):
    """
    get the l2 error based on an EM algorithm

    Parameters
    ----------
        pref : tuple
            a tuple with length n like (0, 1, 3, 2), which means 0 > 1 > 3 > 2
        features : array
            features of each alternative

    Returns
    -------
        error : float
            the estimated l2 error
    """
    inst = ApproxRUM(features)
    inst.EM_algorithm(pref, _mixture(), _max_iter(),fixed_effects)
    beta_est_mat = inst.beta_mat_list[-1]
    print(beta_est_mat)
    lambda_est = inst.lambda_list[-1]
    print(lambda_est)
    print(sum(lambda_est))
    rho_est = inst.get_rho_from_mixed_logit(lambda_est, beta_est_mat,fixed_effects)
    error = inst.l2_distance(inst.rho, rho_est)
    return pref, error





if __name__ == "__main__":
    # Alternatives

    ranking_list = list(permutations((0, 1, 2, 3)))
    # Outcomes depend on the initialization of parameters
    # Compute error for each of the rankings many times to avoid failure of optimization that is caused by bad initial points


    num_init_choice = _n_starts()

    output = []
    # fish data
    features = _load_features()
    fixed_effects = np.zeros(features.shape[0])
#####################################################
######if order=2, use second order polynomial########
#####################################################
    if order==2:
      features=np.concatenate((features,np.square(features)),axis=1)

    if _standardise():
        features = (features - features.mean(0)) / features.std(0)
        print('features standardised', flush=True)
    print('starts', num_init_choice, 'mixture', _mixture(), 'max_iter', _max_iter(),
          'tol', _tol(), flush=True)

    pref = ranking_list[index]
    
    
#####################################################
######EM Algorithm###################################
#####################################################

    for i in range(num_init_choice):
        result= process_EM(pref, features,fixed_effects)
        output.append(result[1])
        
    print(index, output)

#####################################################
######Save Resuls####################################
#####################################################

    file_name = _out_dir() + '/Fish_em_'+ str(order)+'_'+str(index)+ '.txt'
    # Python's min() never updates past a NaN -- x < nan is False -- so a NaN in the
    # first position discards every later start. Take the minimum over the finite starts,
    # and fall back to NaN only when no start returned one.
    finite = [x for x in output if x == x]
    best = min(finite) if finite else float('nan')
    print('starts finite', len(finite), 'of', len(output), 'best', best, flush=True)
    np.savetxt(file_name, np.array([pref[0],pref[1],pref[2],pref[3],best]))