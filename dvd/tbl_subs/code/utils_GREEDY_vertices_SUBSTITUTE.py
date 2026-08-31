import numpy as np
import pandas as pd
from itertools import combinations, permutations
from scipy.optimize import minimize, root
from joblib import Parallel, delayed
from collections import defaultdict
from collections.abc import Iterable
import sys

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
        self.cal_D =[self.X]  + list(combinations(self.X, 2)) + list(combinations(self.X, 3)) 
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
        
    def set_rho_zero(self, target_pref):
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
                self.rho[(x, D)] = 0
        self.rho_set = dict([(D, np.array([self.rho[(x, D)] for x in D])) for D in self.cal_D])
    
    def substitution_difference(self, rho_2, x, x_r):
        """
        Parameters
        ----------
            rho_2 : dict
                a stochastic choice function
            x : int
                choice
            x_r : int
                removed choice
        Returns
        -------
            difference : float
        """
        D_ur = (0, 1, 2, 3)
        l = list(D_ur)
        l.remove(x_r)
        D_r = tuple(l)
        
        val = rho_2[(x, D_ur)]
        for i in D_r:
            if i != x:
                val += rho_2[(i, D_r)]
        return val**2

    def get_rho_from_logit(self, beta, fixed_effects=0):
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
            probs = self.logit(D, beta, fixed_effects[tuple([D])])
            for i, x in enumerate(D):
                rho_est[(x, D)] = probs[i]
        return rho_est

    def greedy_algorithm(self, target_pref, iter_num=100, target_pref_2=None, weight=None, fixed_effects=0, x=1, x_r=0,):
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
            x : int
                choice
            x_r : int
                removed choice

        Returns
        -------
            error: L2 error

        Notes
        -----
            If you want to approximate a mixture point of two vertices, you can specify the other point in target_pref_2.
            In this case, you should also indicate the mixture weight.
        """
        # set the target preference
        self.set_rho_zero(target_pref)

        alpha = lambda n: 1 / (n + 1)

        self.z_list = [dict([(D, np.zeros(len(D))) for D in self.cal_D])]
        self.alpha_list = []
        self.beta_list = []
        self.rho_list = []
        self.l2_error_list = [np.inf]
        self.subs_diff_list = [np.inf]
       
        for n in range(iter_num):
            if n % 100 == 0:
                print(n)
            # z^{n-1}
            rho_est_set_prev = self.z_list[-1]
            rho_est_prev = dict(sum([[((x, D), prob) for x, prob in zip(D, rho_est_set_prev[D])] for D in self.cal_D], []))
            min_val_temp = 100000
           
            for j in range(n+1):
                target_rho_set = dict([(D, (self.rho_set[D] - (1 - alpha(j)) * rho_est_set_prev[D])/alpha(j)) for D in self.cal_D])

            # set target
            # compute beta of z_n
                beta = self.update_beta_substitute(target_rho_set, alpha(j),fixed_effects, x, x_r)
            # z_n
                rho_est_curr = self.get_rho_from_logit(beta, fixed_effects)
                rho_est_set_curr = dict([(D, np.array([rho_est_curr[(x, D)] for x in D])) for D in self.cal_D])

            # z^n
                rho_set_est = dict([(D, (1 - alpha(j)) * rho_est_set_prev[D] + alpha(j) * rho_est_set_curr[D]) for D in self.cal_D])
                rho_est = dict(sum([[((x, D), prob) for x, prob in zip(D, rho_set_est[D])] for D in self.cal_D], []))
                
            # progress just for log
                diff = self.substitution_difference(rho_est, x, x_r)
                if diff <= min_val_temp:
                    min_val_temp = diff
                    rho_set_est_min = dict([(D, (1 - alpha(j)) * rho_est_set_prev[D] + alpha(j) * rho_est_set_curr[D]) for D in self.cal_D])
                    
            # store parameter records
            self.z_list.append(rho_set_est_min)
            self.alpha_list.append(alpha(n))
            self.beta_list.append(beta)
            # print('Difference is')
            # print(diff)
            self.subs_diff_list.append(diff)

        # calculate lambdas
        self.lambda_list = []
        for n in range(iter_num):
            lam = self.alpha_list[n] * np.prod([1 - alpha for alpha in self.alpha_list[n+1:]])
            self.lambda_list.append(lam)

    def update_beta_substitute(self, target_rho_set, n, fixed_effects, x, x_r):
        """
        Parameters
        ----------
            target_rho_set : dict
                a dictionary of which keys are choice sets
            fixed_effects : len(X)-dim array
            x : int
                choice
            x_r : int
                removed choice

        Returns
        -------
            beta_new : K-dim array

        """
        if not isinstance(fixed_effects, Iterable):
            fixed_effects = np.zeros(len(self.X))
            
        def eq(beta, x, x_r):
            D_ur = (0, 1, 2, 3)
            l = list(D_ur)
            l.remove(x_r)
            D_r = tuple(l)
            
            prob_temp_ur=self.logit(D_ur, beta, fixed_effects[tuple([D_ur])])
            prob_temp_r=self.logit(D_r, beta, fixed_effects[tuple([D_r])])
            
            i = D_ur.index(x)
            val = prob_temp_ur[i]
            for i in range(len(D_r)):
                if D_r[i] != x:
                    val += prob_temp_r[i]
            return val**2
        
        beta_init = np.random.normal(scale=1, size=self.K)
        method = "BFGS"
        #method = 'SLSQP'
        res = minimize(lambda beta: eq(beta,x,x_r),beta_init,method=method)
        beta_new = res.x
        return beta_new

def process_greedy(pref, features, fixed_effects=0, x=1, x_r=0):
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
        x : int
            choice
        x_r : int
            removed choice

    Returns
    -------
        error : float
            the estimated l2 error
    """
    inst = ApproxRUM(features)
    inst.greedy_algorithm(pref, 1000, fixed_effects=fixed_effects, x=x, x_r=x_r)
    beta_est = inst.beta_list[-1]
    rho_est = inst.get_rho_from_logit(beta_est)
    diff = (1 - np.sqrt(np.array(inst.subs_diff_list[-2:])))
    return pref, diff
