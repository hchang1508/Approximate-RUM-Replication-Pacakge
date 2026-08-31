"""tbl:in-sample_d1 / tbl:out-sample -- the "our method" row: a 4-mixture mixed logit
fitted by the EM algorithm, no fixed effects, on the single choice set D = {J}.

    python code/approx_rum_EM_fish_inout.py <split_index>   # from the exhibit directory

Reads the aggregated characteristics and choice probabilities that code/in_out_fit.R
wrote for the same split, estimates on the training half, and writes both rows of
predictions and l2 errors to output/raw/em_<split_index>.csv.

The out-of-sample row evaluates the estimated (lambda, beta) at the TEST half's
average characteristics and compares with the test half's choice probabilities,
matching what run.alt.regressions.macro does for every benchmark model.

beta is initialised from np.random.normal under a deterministic seed, so the run
is reproducible from the split index.

M = 4 mixtures, max_iter = 1000, tol 1e-6 on the l2 error, cal_D = [X] (the single
full menu), no fixed effects. With one choice set the class's l2_distance divides
by count = 1, so the reported error is the plain Euclidean ||rho_hat - rho||_2 --
the same metric the R benchmarks use via dist.prob.
"""
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
from itertools import combinations, permutations
from scipy.optimize import minimize, root
from tqdm import tqdm
from joblib import Parallel, delayed
from collections import defaultdict
from collections.abc import Iterable

import os
import warnings
import sys
warnings.simplefilter('ignore')

index1 = int(sys.argv[1])

EM_INPUT_DIR = os.environ.get('EM_INPUT_DIR', 'output/em_input')
OUT_DIR      = os.environ.get('OUT_DIR', 'output/raw')

file_train_data = os.path.join(EM_INPUT_DIR, 'data.train.'   + str(index1) + '.1' + '.csv')
file_train_CP   = os.path.join(EM_INPUT_DIR, 'data.train.CP' + str(index1) + '.1' + '.csv')
file_test_data  = os.path.join(EM_INPUT_DIR, 'data.test.'    + str(index1) + '.1' + '.csv')
file_test_CP    = os.path.join(EM_INPUT_DIR, 'data.test.CP'  + str(index1) + '.1' + '.csv')

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
        self.cal_D = [self.X] # + list(combinations(self.X, 2)) + list(combinations(self.X, 3)) 
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
        gamma_mat = weighted_prob_vec / weighted_prob_vec.sum(axis=0)   # M \times len(D)
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

    def EM_algorithm(self, target_pref, M, max_iter_num, fixed_effects,target_pref_2=None, weight=None):
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
        if target_pref_2 and weight:
            self.set_rho_from_pref_mixture(target_pref, target_pref_2, weight)
        else:
            self.set_rho_from_pref(target_pref)

        print('Target_prob',self.rho_set)
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
            if abs(error - self.l2_error_list[-2]) < 1e-6:
                print(i)
                print('tolerance reached')
                break

            # update parameters
            lambda_old = lambda_new
            beta_old_mat = beta_new_mat

    def set_rho_from_pref(self, target_pref):
        """
        Parameters
        ----------
            target_pref : tuple
                a tuple with length n like (0, 1, 3, 2), which means 0 > 1 > 3 > 2
        """
        self.rho = dict()
       #rank_array = np.argsort(target_pref)
        for D in self.cal_D:
            rho = {}
            #min_rank = rank_array[tuple([D])].min()
            for x in D:
                self.rho[(x, D)] =target_pref[x]
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
            count+=1
            for x in D:
                squered_error_total += (rho_1[(x, D)] - rho_2[(x, D)]) ** 2
        error = (squered_error_total ) ** 0.5/ count
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
        print(fixed_effects)
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


def process_EM(pref, features,fixed_effects,target_pref_2=None, weight=None):
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
    inst.EM_algorithm(pref, 4, 1000,fixed_effects,target_pref_2, weight)
    beta_est_mat = inst.beta_mat_list[-1]
    lambda_est = inst.lambda_list[-1]
    rho_est = inst.get_rho_from_mixed_logit(lambda_est, beta_est_mat,fixed_effects)
    print(rho_est)
    print(inst.rho)
    error = inst.l2_distance(inst.rho, rho_est)
    return pref, error,lambda_est,beta_est_mat





def read_features(path):
    """average characteristics X1.0, X0.1 as written by gen.macro.data + write.csv"""
    # columns are "", alt, chid, X1.0, X0.1; rows in the order beach, boat, charter, pier
    dat = pd.read_csv(path)
    assert list(dat.columns[3:5]) == ['X1.0', 'X0.1'], dat.columns.tolist()
    assert list(dat.iloc[:, 1]) == ['beach', 'boat', 'charter', 'pier'], dat.iloc[:, 1].tolist()
    return dat.iloc[:, [3, 4]].to_numpy()


def read_CP(path):
    """empirical choice probabilities as written by compute.emp.CP + write.csv"""
    cp = pd.read_csv(path).iloc[:, 1].to_numpy()
    assert abs(cp.sum() - 1) < 1e-8, cp
    return cp


if __name__ == "__main__":
    # Alternatives
    # Outcomes depend on the initialization of parameters
    # Compute error for each of the rankings many times to avoid failure of optimization that is caused by bad initial points

    # Deterministic given the split index: beta is initialised at random below.
    np.random.seed(index1)

    # fish data
    features      = read_features(file_train_data)
    features_test = read_features(file_test_data)

    #direct optimizing fixed effects
    #fixed_effects = np.diag([1,1,1,1])
    #features = np.concatenate((features,fixed_effects[:,0:3]),axis=1)

    fixed_effects = np.zeros(4)

    pref      = read_CP(file_train_CP)
    pref_test = read_CP(file_test_CP)

    EM_result = process_EM(pref, features, fixed_effects)
    _, error_in, lambda_est, beta_est_mat = EM_result

    # In-sample predictions: the mixture evaluated at the TRAINING characteristics.
    inst_in  = ApproxRUM(features)
    D        = inst_in.cal_D[0]
    pred_in  = inst_in.mixed_logit(D, beta_est_mat, lambda_est, fixed_effects[list(D)])

    # Out-of-sample predictions: the same estimates evaluated at the TEST
    # characteristics, compared with the test half's choice probabilities. This
    # mirrors run.alt.regressions.macro, which predicts with dat.test.macro and
    # scores against emp.CP.test.
    inst_out  = ApproxRUM(features_test)
    pred_out  = inst_out.mixed_logit(D, beta_est_mat, lambda_est, fixed_effects[list(D)])
    error_out = float(np.sqrt(np.square(pred_out - pref_test).sum()))

    # Guard: process_EM's error must be the l2 distance of pred_in from the training CP.
    assert abs(error_in - np.sqrt(np.square(pred_in - pref).sum())) < 1e-8, error_in

    alts = ['beach', 'boat', 'charter', 'pier']
    out  = pd.DataFrame([
        dict(split=index1, spec='our_method', sample='in',
             **dict(zip(alts, pred_in)),  loss=error_in),
        dict(split=index1, spec='our_method', sample='out',
             **dict(zip(alts, pred_out)), loss=error_out),
    ])

    os.makedirs(OUT_DIR, exist_ok=True)
    out.to_csv(os.path.join(OUT_DIR, 'em_%d.csv' % index1), index=False)

    # Parameter estimates, for inspection; not consumed by collect.R.
    np.savetxt(os.path.join(OUT_DIR, 'em_lambda_%d.csv' % index1), np.array(lambda_est), delimiter=',')
    np.savetxt(os.path.join(OUT_DIR, 'em_beta_%d.csv'   % index1), np.array(beta_est_mat), delimiter=',')

    print(out.to_string(index=False))

 
