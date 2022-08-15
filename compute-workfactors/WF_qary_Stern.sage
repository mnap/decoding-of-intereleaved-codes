"""
Compute WF of q-ary Stern.

REFERENCE PAPER: 
Peters, Christiane. 2010. ‘Information-Set Decoding for Linear Codes over Fq’. In Post-Quantum Cryptography, edited by Nicolas Sendrier, 81–94. Lecture Notes in Computer Science. Berlin, Heidelberg: Springer. https://doi.org/10.1007/978-3-642-12929-2_7.

The paper presents the q-ary Stern algorithm and also various improvements of this algorithm.
Depending on which improvements we consider, the computation of the WF can become very expensive.

It seems the paper presents the following improvements:
 - reusing additions (sec 4)
 - collisions (sec 4)
 - reusing parts of information sets and precomputations (sec 4)
 - overlapping sets (sec 7)

WF_qary_Stern
 - This function computes the WF for the q-ary Stern algorithm where only the first two improvements
listed above are included.
 - This algorithm has two internal parameters p_int and l_int. While p_int is optimized exhaustively,
a heuristic is used for l_int.
 - Since k can be odd, k/2 is replaced by k//2 (i.e. floor(k/2)). However this isn't exactly correct.
For example if k is odd the correct number of expected collisions is proportional to:
(k//2 choose p)*((k//2)+1 choose p) insead of (k//2 choose p)^2
"""

import logging
logging.basicConfig(level=logging.DEBUG)
debug = logging.debug
logging.disable(level=logging.DEBUG) # comment out to enable debug messages


def WF_qary_Stern(q, n, k, t):
    best_WF_bits = float('inf')
    best_param = dict()
    # internal parameters: p_int, l_int
    for p_int in range(0, t//2 + 1): # loop over p_int
        l_int_guess = round(log(binomial(k//2, p_int), q) + p_int*log(q - 1, q))
        l_int_MIN = max(0, floor(l_int_guess - 0.25*t))
        l_int_MAX = min(n-k, ceil(l_int_guess + 0.25*t))
        for l_int in range(l_int_MIN, l_int_MAX + 1):
            WF_bits, exp_iterations_bits = WF_qary_Stern__helper(q=q, n=n, k=k, t=t, p_int=p_int, l_int=l_int)
            if best_WF_bits > WF_bits:
                best_WF_bits = WF_bits
                best_param['p_int'] = p_int
                best_param['l_int'] = l_int
                best_param['l_int_guess'] = l_int_guess # for reference
                best_param['exp_iterations_bits'] = exp_iterations_bits # for reference

    return best_WF_bits, best_param


def WF_qary_Stern__helper(q, n, k, t, p_int, l_int):
    one_iteration_cost = (n - k)^2*(n + k) \
                        + ((k//2 - p_int + 1) + 2*binomial(k//2, p_int)*(q - 1)^p_int)*l_int \
                        + (q/(q - 1))*(t - 2*p_int + 1)*2*p_int *(1 + (q - 2)/(q - 1)) \
                        * (binomial(k//2, p_int)^2*(q-1)^(2*p_int))/q^l_int
    success_probability = binomial(k//2, p_int)^2 * binomial(n - k - l_int, t - 2*p_int)/binomial(n, t)
    exp_iterations_bits = round(log(1/success_probability, 2), ndigits=2)
    WF = N(one_iteration_cost/success_probability)
    WF_bits = log(WF, 2)
    debug(f'{q=}, {n=}, {k=}, {t=}, {p_int=}, {l_int=} --> {WF_bits=}')
    return WF_bits, exp_iterations_bits
