from sage_import import sage_import
sage_import('WF_qary_Stern_Dumer', fromlist=['WF_qary_Stern_Dumer'])
sage_import('interleaved_workfactors', fromlist=['get_random_prange_wf', 'get_int_prange_wf'])


def log2(x):
    return log(x, 2)


def test__yield_params(count=None):
    params = [
        "\n\n-------------------WF_DESIGN=128------------\n",
        dict(n=2130, r=100, m=8, k=1330, q=3, ell=7, SL=128, d_E= 59, t_pub=131),
        dict(n=1580, r= 90, m=6, k=1040, q=4, ell=7, SL=128, d_E= 68, t_pub=105),
        dict(n=1290, r=100, m=5, k= 790, q=5, ell=7, SL=128, d_E= 73, t_pub=109),
        "\n\n-------------------WF_DESIGN=256------------\n",
        dict(n=4300, r=180, m=8, k=2860, q=3, ell=7, SL=256, d_E=123, t_pub=236),
        dict(n=3760, r=240, m=7, k=2080, q=4, ell=7, SL=256, d_E=135, t_pub=280),
        dict(n=3200, r=200, m=6, k=2000, q=5, ell=7, SL=256, d_E=121, t_pub=218),
    ]

    for param in params[:count]:
        yield param


if __name__ == '__main__':
    for prm in test__yield_params():
        if isinstance(prm, str):
            print(prm)
            continue
        q, ell, n, m, r, t, SL = [prm[x] for x in ['q', 'ell', 'n', 'm', 'r', 't_pub', 'SL']]
        k = n - m*r
        assert k == prm['k']
        d_E = prm['d_E']

        try:
            # we compute the WF of Stern using the Stern-Dumer variant
            int_stern_dumer_wf_log, int_stern_dumer__best_prm = 0, 0
            int_stern_dumer_wf_log, int_stern_dumer__best_prm = WF_qary_Stern_Dumer(q=q, n=n, k=k+ell, t=d_E) # in bits
        except Exception:
            traceback.print_exc()

        random_prange_wf = get_random_prange_wf(q=q, n=n, k=k, t=t)
        int_prange_wf = get_int_prange_wf(q=q, n=n, k=k, t=t, ell=ell)
        int_stern_dumer_wf = 2^int_stern_dumer_wf_log

        ########
        # OUTPUT
        ########
        print("q={q}, ell={ell}, n={n}, k={k}, t={t},"
              .format(q=q, ell=ell, n=n, k=k, t=t))
        print("int_stern_dumer__best_prm", int_stern_dumer__best_prm)

        header_column = [
            " ",
            "CF-based Stern-Dumer WF",
            "IntPrange WF",
            "Random Prange WF",
        ]
        results = [
            int_stern_dumer_wf*log2(q),
            int_prange_wf*log2(q),
            random_prange_wf*log2(q),
        ]

        assert len(results) + 1 == len(header_column)
        results = [N(x) for x in results]
        results_log2 = [log2(x) for x in results]
        results = ["{:.2e}".format(N(x)) for x in results]
        results_log2 = ["{:6.2f}".format(N(x)) for x in results_log2]

        pretty_output = table(columns=[["log2(X)"] + results_log2,
                                       ["X"] + results],
                              header_row=True,
                              header_column=header_column,
                              frame=False)
        print(pretty_output)
        print()
