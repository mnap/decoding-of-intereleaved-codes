"""
Compute WF of q-ary Stern-Dumer.

Heuristics:
 - it seems the WF is convex in l_int, so one can take advantage of that to find the minimum WF over l_int
 - p/2 is replaced by by p//2, so we search over only even p
 in WF_qary_Stern_Dumer__helper:
   - we do (k+l)/2 --> (k+l)//2 and p/2 --> p//2
   - we include the cost for Gaussian elimination (it is needed at least for the case p_int=0 so that we get a non-zero WF)
"""

import logging
logging.basicConfig(level=logging.DEBUG)
debug = logging.debug
logging.disable(level=logging.DEBUG) # comment out to enable debug messages


def WF_qary_Stern_Dumer(q, n, k, t):
    best_WF_bits = float('inf')
    best_param = dict()

    # internal parameters: p_int, l_int
    # where 0 <= p_int <= t and 1 <= l_int <= n-k-(t-p_int) # TODO verify
    for p_int in range(0, t+1, 2): # loop over p_int
        debug("=========new p_int==========")
        MIN_L_INT = 1
        MAX_L_INT = n-k-t+p_int
        # for reference we will record the range of values that we tried for l_int
        max_l_int_tried = float('-inf')
        min_l_int_tried = float('inf')
        WF_bits_at_max_l_int_tried = float('inf')
        WF_bits_at_min_l_int_tried = float('inf')
        # using the formula from Peters' paper
        # note we have p_int//2 instead of p_int because in Peter's algorithm, p_int corresponds to
        # half of the p_int of this algorithm
        l_int_guess = round(log(binomial(k//2, p_int//2), q) + (p_int//2)*log(q - 1, q))
        l_int_guess = max(MIN_L_INT, l_int_guess)
        l_int_guess = min(MAX_L_INT, l_int_guess)
        # search optimal value of l_int assuming the WF is convex in l_int
        # method:
        # starting at l_int=l_int_guess, increment l_int until WF starts increasing
        # then starting at l_int=l_int_guess-1, decrement l_int until WF starts increasing
        # finally take the minimum of all the WF's we just computed
        temp_best_WF_bits = float('inf')
        temp_best_param = dict()
        for mode in ["increasing", "decreasing"]:
            if mode == "increasing":
                l_int = l_int_guess
                var_change = 1
            else:
                l_int = l_int_guess - 1
                var_change = -1
            while True:
                if not (MIN_L_INT <= l_int and l_int <= MAX_L_INT):
                    break
                WF_bits = WF_qary_Stern_Dumer__helper(q=q, n=n, k=k, t=t, p_int=p_int, l_int=l_int)
                if l_int > max_l_int_tried:
                    max_l_int_tried = l_int
                    WF_bits_at_max_l_int_tried = WF_bits
                if l_int < min_l_int_tried:
                    min_l_int_tried = l_int
                    WF_bits_at_max_l_int_tried = WF_bits
                if temp_best_WF_bits > WF_bits:
                    temp_best_WF_bits = WF_bits
                    # l_int_guess is included just for our reference
                    temp_best_param = dict(p_int=p_int, l_int=l_int, l_int_guess=l_int_guess)
                    l_int += var_change
                else:
                    break
        debug("tried l_int in [{0}, {1}]".format(min_l_int_tried, max_l_int_tried))

        if best_WF_bits > temp_best_WF_bits:
            best_WF_bits = temp_best_WF_bits
            best_param = temp_best_param

    return best_WF_bits, best_param


def WF_qary_Stern_Dumer__helper(q, n, k, t, p_int, l_int):
    w = t

    num = binomial(n, w)
    den = binomial((k+l_int)//2, p_int//2)^2 * binomial(n-k-l_int, w-p_int)
    term_1 = num/den

    term_21 = 2*l_int*intermediate_sum(q, (k+l_int)//2, p_int//2)
    term_22 = binomial((k+l_int)//2, p_int//2)^2 * (q-1)^p_int/q^l_int * q/(q-1) * (w-p_int+1)*p_int
    #debug("{:.2f} {:.2f} {:.2f}".format(float(term_1), float(term_21), float(term_22)))

    term_0 = n*(n-k-l_int)*(n-k) # cost for Gaussian elimination

    WF = term_1*(term_0 + term_21 + term_22)
    WF_bits = N(log(WF, 2))
    debug(f'{q=}, {n=}, {k=}, {t=}, {p_int=}, {l_int=} --> {WF_bits=}')
    assert WF_bits != -Infinity # for debugging
    return WF_bits


def intermediate_sum(q, b, s):
    total = 0
    for i in range(1, s+1):
        total += binomial(b, i)*(q-1)^i
    return total
