import os
import traceback # for debugging
from pprint import pprint
from functools import lru_cache
from sage.misc.table import table
from sage_import import sage_import


def get_gaussian_elimination_cost(n):
    """Return cost of reducing a size (n x n) matrix A to row echelon form.
    This also equals the cost of solving Ax = b for x (ignoring smaller order terms).
    """
    return n^3
    

def prob_lin_ind__helper(q, nrows, ncols, jstart=0):
    '''Return prod (1 - q^(j - ncols)) from j = jstart to j = nrows - 1.

    If jstart = 0, this equals the probability that a random matrix in F_q^(nrows x ncols) is
    full-rank.
    '''
    prod = 1
    for j in range(jstart, nrows):
        prod *= 1 - q^(j-ncols)
    return prod


def get_prange_Psucc(n, k, t):
    """Return Prange's success probability."""
    return binomial(n-t, k)/binomial(n, k)


def get_random_prange_Psucc(q, n, k, t):
    """Return Random Prange's success probability.

    That is, pick a random nonzero codeword in <R> and run one iteration of Prange on it.
    (As usual, the error matrix is assumed to be chosen at random.)"
    """
    P_succ = 0
    for w in range(0, t+1):
        first_term = binomial(t, w)*((q-1)^w)/(q^t)
        second_term = binomial(n-w, k)/binomial(n, k)
        P_succ += first_term*second_term
    return P_succ


def get_random_prange_wf(q, n, k, t):
    """Return WF of Random Prange."""
    return get_gaussian_elimination_cost(k)/get_random_prange_Psucc(q=q, n=n, k=k, t=t)


def prob_GprimeJ_is_rank_deficient_given_etc(q, ell, k, p_hat):
    '''Return probability that (G_prime)_J has a rank deficiency of p_hat given that G_J is
    full-rank and E_J has linearly independent rows.

    (This is in the context of intPrange.)
    '''
    term1 = 1
    for j in range(0, k):
        term1 *= (q^k - q^j)/(q^(k+ell) - q^j)
    term2 = gaussian_binomial(n=ell, k=p_hat, q=q)
    term3 = q^((ell - p_hat)*(k - p_hat))
    term4 = gaussian_binomial(n=k, k=k - p_hat, q=q)

    return N(term1*term2*term3*term4)


def get_int_prange_Pfail(q, n, k, t, ell):
    """Return probability that an iteration of IntPrange leads to a failure (i.e. the rank-check
    step succeeds succeeds but the error vector is not found from the last step).
    In other words, this computes the probability that
       1. G'_J is rank deficient and
       2. E_J has linearly independent rows and
       3. G_J is full rank (i.e. J is an information set)
    is true.

    Assumes random error matrices."""
    Psucc = get_int_prange_Psucc(q=q, n=n, k=k, t=t, ell=ell)
    P_GJ_has_full_rank = prob_lin_ind__helper(q=q, nrows=k, ncols=k+ell)
    main_term = P_GJ_has_full_rank - Psucc
    next_term = 1 - (prob_lin_ind__helper(q=q, nrows=k+ell, ncols=k+ell, jstart=k) \
                     /prob_lin_ind__helper(q=q, nrows=ell, ncols=k+ell, jstart=0))
    ret = N(next_term*main_term)
    return ret


@lru_cache
def get_int_prange_Psucc(q, n, k, t, ell):
    """IntPrange's success probability when the error matrix is chosen at random.

    In other words, computes the probability that
       1. G'_J is rank deficient and
       2. E_J has linearly dependent rows and
       3. G_J is full rank (i.e. J is an information set)
    is true.
    """
    P_GJ_has_full_rank = prob_lin_ind__helper(q=q, nrows=k, ncols=k+ell)
    P_succ = 0
    for i in range(0, t+1):
        first_term = binomial(n-t, k+ell-i)*binomial(t, i)/binomial(n, k+ell)
        second_term = 1 - prob_lin_ind__helper(q=q, nrows=ell, ncols=i)
        P_succ += first_term*second_term
    P_succ *= P_GJ_has_full_rank
    return P_succ


def get_int_prange_last_step_wf(q, k, p_hat):
    """Retrun cost of performing last step of intPrange.

    In the last step we do q^p_hat checks where p_hat is the rank deficiency of G'_J."""
    return (q^p_hat) * get_gaussian_elimination_cost(k)


def get_int_prange_wf(q, n, k, t, ell):
    '''Return WF of intPrange.'''
    # mainterm1
    WF_rank_check = get_gaussian_elimination_cost(k+ell)
    P_succ = get_int_prange_Psucc(q=q, n=n, k=k, t=t, ell=ell)
    P_fail = get_int_prange_Pfail(q=q, n=n, k=k, t=t, ell=ell)
    mainterm1 = WF_rank_check*(1 - P_fail - P_succ)/(1 - P_succ)

    # mainterm2
    P_GJ_has_full_rank = prob_lin_ind__helper(q=q, nrows=k, ncols=k+ell)
    term1 = (1 - P_succ/P_GJ_has_full_rank)/(1 - P_succ)
    term2 = 1
    for j in range(0, k):
        term2 *= (q^k - q^j)/(q^(k + ell))
    term3 = 0
    for p_hat in range(1, ell + 1):
        term31 = gaussian_binomial(n=ell, k=p_hat, q=q)
        term32 = q^((ell - p_hat)*(k - p_hat))
        term33 = gaussian_binomial(n=k, k=k - p_hat, q=q)
        term34 = WF_rank_check + get_int_prange_last_step_wf(q=q, k=k, p_hat=p_hat)
        term3 += term31*term32*term33*term34
    mainterm2 = term1*term2*term3

    ret = (mainterm1 + mainterm2)/P_succ
    return N(ret)
