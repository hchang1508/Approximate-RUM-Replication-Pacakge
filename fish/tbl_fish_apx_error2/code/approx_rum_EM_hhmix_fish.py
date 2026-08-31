import numpy as np
import pandas as pd
from scipy.optimize import minimize, root
from collections import defaultdict
from collections.abc import Iterable
from itertools import combinations, permutations

import warnings
import sys
import random
import os
def _load_features(default='../../data/fish/fish_J4_K2.csv'):
    path = sys.argv[4] if len(sys.argv) > 4 else os.environ.get('FISH_DATA', default)
    df = pd.read_csv(path)
    print('characteristics from', path, flush=True)
    print(df.to_string(index=False), flush=True)
    return df.drop(columns=['alt_id']).to_numpy(dtype=float)

def _out_dir():
    d = os.environ.get('OUT_DIR', 'output/raw')
    os.makedirs(d, exist_ok=True)
    return d

def _n_starts(default=10):
    # Random initialisations of EM per grid point. The seed is
    # drawn once per (case, degree, chunk) and consumed in sequence, so raising this leaves
    # the first 10 starts identical and can only lower the reported minimum.
    return int(os.environ.get('EM_STARTS', default))

def _tol(default=1e-4):
    # Termination threshold on the absolute change in error between consecutive EM
    # iterations. This is an absolute test, so on a badly-scaled design it fires while the
    # fit is still far from optimal; submit_em.sh tightens it to 1e-8.
    return float(os.environ.get('EM_TOL', default))

def _standardise():
    # Centre and scale the feature columns. This is a reparametrisation: utilities shift by
    # a constant within each choice set (which cancels in the logit) and beta rescales, so
    # the attainable choice probabilities -- and hence the approximation error -- are
    # unchanged. It only conditions the iteration path. Off by default.
    return os.environ.get('EM_STANDARDIZE', '0') not in ('0', '', 'false', 'False')

def _chunk(default=300):
    # Grid points per array task. Smaller chunks give shorter tasks, which backfill better
    # and lose less to a preemption or a node failure; submit_em.sh uses 25.
    return int(os.environ.get('EM_CHUNK', default))
warnings.simplefilter('ignore')

# Generator behind every EM initial point. It is reseeded from the run's derived seed in
# __main__; the value here only covers imports that never reach that block.
#
# Until 2026-08-16 this was a hardcoded default_rng(0) and nothing reseeded it. Seeding at
# the top of __main__ used np.random.seed(), which drives the legacy global RNG and has no
# effect on a Generator object, so every run drew the same initial points regardless of
# SEED_BASE. Restarts taken *within* a task were unaffected -- the generator is consumed in
# sequence, so successive starts differ -- but restarts taken as separate replicate
# processes collapsed onto one another. Set EM_RNG_SEED=0 to reproduce that behaviour.
rng = np.random.default_rng(0)

def _seed_rng(seed):
    """Reseed the module-level generator used for EM initial points."""
    global rng
    rng = np.random.default_rng(seed)
    return rng

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
        self.cal_D = [self.X]  + list(combinations(self.X, 2)) + list(combinations(self.X, 3)) 
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

        #print('Target_prob',self.rho_set)
        # initialize parameters
        lambda_old = np.zeros(M) + 1 / M
        beta_old_mat = rng.normal(scale=0.5, size=(M, self.K))
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
                break
            
            

            # update parameters
            lambda_old = lambda_new
            beta_old_mat = beta_new_mat

        if i == max_iter_num - 1:
            # range(max_iter_num) leaves i at max_iter_num-1, so the old `i==max_iter_num`
            # test never fired and a capped-out fit was indistinguishable from a converged
            # one in the logs.
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
        print('fixed_effects')
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

        if not isinstance(fixed_effects, Iterable):
            fixed_effects = np.zeros(len(self.X))

        rho_est = dict()
        for D in self.cal_D:

            probs = self.mixed_logit(D, beta_mat, lambda_vec, fixed_effects[list(D)])
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

        beta_init = rng.normal(scale=0.1, size=self.K)
        method = "BFGS"
        res = minimize(eq, beta_init, bounds=[(-20, 20) for _ in range(self.K)], method=method)
        beta_new = res.x
        return beta_new


def process_EM(pref, features,fixed_effects,target_pref_2=None, num_mixture=4,weight=None):
 
    
    #Initiate an ApproxRUM instance
    inst = ApproxRUM(features)

    #EM algorithm
    inst.EM_algorithm(pref, num_mixture, 1000,fixed_effects,target_pref_2, weight)
    
    #estimated coefficients
    beta_est_mat = inst.beta_mat_list[-1]
 
    #estimated mixture weights
    lambda_est = inst.lambda_list[-1]
    
    #approximated choice probabilities
    rho_est = inst.get_rho_from_mixed_logit(lambda_est, beta_est_mat,fixed_effects)
    
    #errors measure in the l2/|choice sets| distance
    error = inst.l2_distance(inst.rho, rho_est)
    
    return pref, error,lambda_est,beta_est_mat


if __name__ == "__main__":


    features = _load_features()
    # a grid of fixed effects
    fixed_effects_dict = [[x,y,z] for x in range(-10,11) for y in range(-10,11) for z in range(-10,11) ]

 
    index = int(sys.argv[1]) #starting point of the fixed effects
    case  = int(sys.argv[2]) #which half-half mixture to be computed
    d = int(sys.argv[3]) #degree of the polynomial

    _SEED = int(os.environ.get('SEED_BASE', 20240808)) + 100000 * d + 1000 * case + index
    np.random.seed(_SEED)
    random.seed(_SEED)
    # EM initial points come from the module-level Generator, which np.random.seed does not
    # touch. Seed it explicitly, or two runs differing only in SEED_BASE return identical
    # numbers. EM_RNG_SEED=0 restores the pre-2026-08-16 behaviour.
    _RNG_SEED = int(os.environ.get('EM_RNG_SEED', _SEED))
    _seed_rng(_RNG_SEED)
    print('seed', _SEED, 'rng seed', _RNG_SEED, flush=True)
    n_starts = _n_starts()
    print('random starts per grid point', n_starts, flush=True)

   
    #the program will compute for fixed effects indexed from index_start to index_end
    chunk = _chunk()
    index_start = 0 + (index-1)*chunk
    index_end = min(0 + (index)*chunk,len(fixed_effects_dict))

    if d==2:
        features=np.concatenate((features,features**2),1)

    if _standardise():
        features = (features - features.mean(0)) / features.std(0)
        print('features standardised', flush=True)
    print('chunk', chunk, 'range', index_start, index_end, 'tol', _tol(), flush=True)

    #for mixture
    if case==1:
        pref_1, pref_2, weight = (0, 1, 2, 3), (3, 2, 1, 0), 0.5 #1234,4321
    elif case==2:
        pref_1, pref_2, weight = (0, 1, 3, 2), (2, 3, 1, 0), 0.5 #1243,3421
    elif case==3:
        pref_1, pref_2, weight = (0, 2, 1, 3), (3, 1, 2, 0), 0.5 #1324,4231
    elif case==4:
        pref_1, pref_2, weight = (0, 3, 1, 2), (2, 1, 3, 0), 0.5 #1423,3241
    elif case==5:
        pref_1, pref_2, weight = (1, 0, 2, 3), (3, 2, 0, 1), 0.5 #2134,4312
    elif case==6:
        pref_1, pref_2, weight = (1, 0, 3, 2), (2, 3, 0, 1), 0.5 #2143,3412

    error=list()
    
    
    #21^3 is the size of the grid
    for i in range(index_start,index_end):

        # fixed effects for this iteration
        fixed_effects= np.array(fixed_effects_dict[i])

        fixed_effects = np.append(fixed_effects,0)
        fixed_effects = np.array(fixed_effects) 
        
        print(i,' Processing:',fixed_effects)
        #number of random initilization
        for j in range(n_starts):
  
            #We set M=18 according to the first paragraph of Section 5
            [_,error_temp,_,_]=process_EM( pref_1,features,fixed_effects,pref_2,num_mixture=18,weight=weight)


            error.append(error_temp)
    #save output            
    error=pd.DataFrame(error)
    output_name=_out_dir() + '/hh_Fish_'+str(case)+'d'+str(d)+'_EM_error_'+str(index_start)+'_'+str(index_end)+'.csv'
    error.to_csv(output_name)

