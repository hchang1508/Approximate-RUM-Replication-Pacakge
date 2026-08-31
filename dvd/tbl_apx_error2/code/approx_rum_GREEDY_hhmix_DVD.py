import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
from itertools import combinations, permutations
from scipy.optimize import minimize, root
from tqdm import tqdm
from joblib import Parallel, delayed
from collections import defaultdict
from collections.abc import Iterable
import sys
import warnings
import time
import os
def _load_features(default='../../data/dvd/dvd_J4_K2.csv'):
    path = sys.argv[4] if len(sys.argv) > 4 else os.environ.get('DVD_DATA', default)
    df = pd.read_csv(path)
    print('characteristics from', path, flush=True)
    print(df.to_string(index=False), flush=True)
    return df.drop(columns=['alt_id']).to_numpy(dtype=float)

def _out_dir():
    d = os.environ.get('OUT_DIR', 'output/raw')
    os.makedirs(d, exist_ok=True)
    return d


import random

warnings.simplefilter('ignore')



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
        self.cal_D = list(combinations(self.X, 2)) + [self.X]  + list(combinations(self.X, 3))
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
            count +=1
            for x in D:
                
                squered_error_total += (rho_1[(x, D)] - rho_2[(x, D)]) ** 2
        #print(count)
        error = (squered_error_total) ** 0.5/count
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

        if not isinstance(fixed_effects, Iterable):
            fixed_effects = np.zeros(len(self.X))
        rho_est = dict()
        for D in self.cal_D:
            probs = self.logit(D, beta, fixed_effects[list(D)])
            for i, x in enumerate(D):

                rho_est[(x, D)] = probs[i]
        return rho_est

    def get_rho_from_mixed_logit(self, lambda_vec, beta_mat, fixed_effects=0):
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

    def greedy_algorithm(self, target_pref, iter_num=100, target_pref_2=None, weight=None, fixed_effects=0):
        """
        Haoge's greedy algorithm

        Parameters
        ----------
            target_pref : tuple
                a tuple with length n like (0, 1, 3, 2), which means 0 > 1 > 3 > 2
            iter_num : int
                number of iterations
            target_pref_2 : tuple, default None
                a tuple with length n like (0, 1, 3, 2), which means 0 > 1 > 3 > 2
            weight : float, default None
                weight for target_pref
            fixed_effects : len(X)-dim array, default 0
                fixed effects

        Returns
        -------
            error: L2 error

        Notes
        -----
            If you want to approximate a mixture point of two vertices, you can specify the other point in target_pref_2.
            In this case, you should also indicate the mixture weight.
        """
        # set the target preference
        if target_pref_2 and weight:
            self.set_rho_from_pref_mixture(target_pref, target_pref_2, weight)
        else:
            self.set_rho_from_pref(target_pref)

        alpha = lambda n: 1 / (n + 1)

        self.z_list = [dict([(D, np.zeros(len(D))) for D in self.cal_D])]
        self.alpha_list = []
        self.beta_list = []
        self.rho_list = []
        self.l2_error_list = [np.inf]
       
        for n in range(iter_num):
            
            print(n)
            # z^{n-1}
            rho_est_set_prev = self.z_list[-1]
            rho_est_prev = dict(sum([[((x, D), prob) for x, prob in zip(D, rho_est_set_prev[D])] for D in self.cal_D], []))
            min_val_temp = 100000
           
            for j in range(n+1):
            # set target
                target_rho_set = dict([(D, (self.rho_set[D] - (1 - alpha(j)) * rho_est_set_prev[D])/alpha(j)) for D in self.cal_D])


            # compute beta of z_n
                beta = self.update_beta(target_rho_set, alpha(j),fixed_effects)

            # z_n
                rho_est_curr = self.get_rho_from_logit(beta, fixed_effects)
                rho_est_set_curr = dict([(D, np.array([rho_est_curr[(x, D)] for x in D])) for D in self.cal_D])

            # z^n
                rho_set_est = dict([(D, (1 - alpha(j)) * rho_est_set_prev[D] + alpha(j) * rho_est_set_curr[D]) for D in self.cal_D])
                rho_est = dict(sum([[((x, D), prob) for x, prob in zip(D, rho_set_est[D])] for D in self.cal_D], []))

            # progress just for log
                error = self.l2_distance(self.rho, rho_est)
               
                if error <= min_val_temp:
                    #print(n,j)
                    #print(error)
                    min_val_temp = error
                    rho_set_est_min = dict([(D, (1 - alpha(j)) * rho_est_set_prev[D] + alpha(j) * rho_est_set_curr[D]) for D in self.cal_D])
                    
            # store parameter records
            self.z_list.append(rho_set_est_min)
            self.alpha_list.append(alpha(n))
            self.beta_list.append(beta)
            self.l2_error_list.append(error)

        # calculate lambdas
        self.lambda_list = []
        for n in range(iter_num):
            lam = self.alpha_list[n] * np.prod([1 - alpha for alpha in self.alpha_list[n+1:]])
            self.lambda_list.append(lam)

    def update_beta(self, target_rho_set, n,fixed_effects=0):
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
        #print(target_rho_set)
        #print(n)
        def eq(beta):
            return sum([np.square(target_rho_set[D] - self.logit(D, beta, fixed_effects[[D]])).sum() for D in self.cal_D])

        beta_init = np.random.normal(scale=1, size=self.K)
        method = "BFGS"
        #method = 'SLSQP'
        res = minimize(eq, beta_init, bounds=[(-20, 20) for _ in range(self.K)], method=method)
        #print(res)
        beta_new = res.x
        return beta_new



def process_greedy(pref, features, fixed_effects=0):
    """
    get the l2 error based on a greedy algorithm

    Parameters
    ----------
        pref : tuple
            a tuple with length n like (0, 1, 3, 2), which means 0 > 1 > 3 > 2
        features : array
            features of each alternative
        fixed_effects : len(X)-dim array, default 0
            fixed effects

    Returns
    -------
        error : float
            the estimated l2 error
    """
    # pref = (2, 0, 1, 3)
    inst = ApproxRUM(features)
    inst.greedy_algorithm(pref, 1000, fixed_effects=fixed_effects)
    beta_est = inst.beta_list[-1]
    rho_est = inst.get_rho_from_logit(beta_est)
    error = inst.l2_error_list[-1]
    return pref, error



if __name__ == "__main__":


    features = _load_features()
        
    # a grid of fixed effects
    fixed_effects_dict = [[x,y,z] for x in range(-10,11) for y in range(-10,11) for z in range(-10,11) ]


    index = int(sys.argv[1])     #index for the fixed effect
    case  = int(sys.argv[2])     #which half-half mixture to be computed
    d = int(sys.argv[3])         #degree of the polynomial
    _SEED = int(os.environ.get('SEED_BASE', 20240808)) + 100000 * d + 1000 * case + index
    np.random.seed(_SEED)
    random.seed(_SEED)
    print('seed', _SEED, flush=True)
    
    
    #if d=2, use second degree polynomial
    if d==2:
        features=np.concatenate((features,features**2),1)
        
    #read to do file
    to_do=pd.read_csv(os.environ.get('TODO_DIR', 'to_do_list_GREEDY') + '/' + str(case)+'_d'+str(d)+'_case_to_do.csv',header=None).to_numpy()
    index = int(to_do[index])
    
   # fixed effects for this iteration
    fixed_effects= np.array(fixed_effects_dict[index])

    fixed_effects = np.append(fixed_effects,0)
    fixed_effects = np.array(fixed_effects) 
    
    #for mixture
    if case==1:
        pref_1, pref_2, weight = (0, 2, 3, 1), (1, 3, 2, 0), 0.5 #1342,2431 
    elif case==2:
        pref_1, pref_2,  weight = (0, 2, 1, 3), (3, 1, 2, 0), 0.5 #1324,4231
    elif case==3:
        pref_1, pref_2, weight = (0, 3, 2, 1), (1, 2, 3, 0), 0.5 #1432, 2341
    elif case==4:
        pref_1, pref_2,  weight = (1, 3, 0, 2), (2, 0, 3, 1), 0.5 #2413,3142
    elif case==5:
        pref_1, pref_2, weight = (1, 0, 3, 2), (2, 3, 0, 1), 0.5 #2143,3412
    elif case==6:
        pref_1, pref_2,  weight = (2, 0, 1, 3), (3, 1, 0, 2), 0.5 #3124,4213


    inst = ApproxRUM(features)
    inst.greedy_algorithm(pref_1, 1000, pref_2, weight, fixed_effects)

    #save output
    output = np.append(fixed_effects,inst.l2_error_list[-1])


    file_name = _out_dir() + '/DVD_hh_gd_' + str(case)+'_'+str(d)+'_'+str(index)+ '.txt'
    np.savetxt(file_name,output,delimiter=",")


