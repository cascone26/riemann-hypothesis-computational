#!/usr/bin/env python3
"""§339a: is reach-to-edge BAND-INDEPENDENT, and does the required exact set grow to the
Riemann-Siegel scale n* ~ gamma/2pi at higher bands (the §337 prediction)?

Bands: [14,44] (RS scale ~3-5), [50,80] (~9-12), [100,130] (~16-21).
Chains: gap(k=3) with exact set = top-R prime powers, R in {3, adaptive}.
All chain quantities are zero-free (analytic side); true margin = lambda_min(pole+arch-prime),
construction validated in §335 (residual 1e-10 vs zero side on [14,44] and [50,80]).
Nyquist edge estimated by Riemann-von Mangoldt: N(T) = (T/2pi)log(T/2pie) + 7/8.
"""
import numpy as np
import l2_edge_sweep as L
from scipy.special import digamma

w = 3.0
# extend the grid for high bands
L.RGRID = np.arange(-170.0, 170.0, 0.002)
L.OMEGA = np.real(digamma(0.25 + 0.5j * L.RGRID)) - np.log(np.pi)

def rvm(T):
    return T / (2 * np.pi) * np.log(T / (2 * np.pi * np.e)) + 7.0 / 8.0

def prime_terms(cs):
    out = {}
    for i, n in enumerate(L.NS):
        if L.LAMW[i] == 0.0:
            continue
        u = L.LOGN[i]
        M = np.zeros((len(cs), len(cs)))
        for j in range(len(cs)):
            for k in range(j, len(cs)):
                v = 2.0 * L.LAMW[i] * L.g_pair(cs[j], cs[k], w, u)
                M[j, k] = M[k, j] = v
        if np.max(np.abs(M)) > 1e-14:
            out[int(n)] = M
    return out

def gap_bound(M, Q, k=3):
    n = M.shape[0]
    k = min(k, n - 1)
    lM, UM = np.linalg.eigh(M)
    Ub, Ur = UM[:, :k], UM[:, k:]
    top = float(np.linalg.eigvalsh(np.diag(lM[:k]) - Ub.T @ Q @ Ub)[0])
    q = float(np.linalg.norm(Ub.T @ Q @ Ur, 2))
    pr = float(np.linalg.eigvalsh(Ur.T @ Q @ Ur)[-1])
    return float(np.linalg.eigvalsh(np.array([[top, -q], [-q, lM[k] - pr]]))[0])

BANDS = [(14.0, 44.0), (50.0, 80.0), (100.0, 130.0)]

for band in BANDS:
    c0, c1 = band
    gmid = (c0 + c1) / 2
    n_eff_edge = None
    print(f"\n=== band [{c0:.0f},{c1:.0f}]  (RS scale gamma/2pi ~ {c0/(2*np.pi):.1f}-{c1/(2*np.pi):.1f}) ===")
    # Nyquist spacing estimate: band width / RvM count in effective band
    Neff = rvm(c1 + 2 * w) - rvm(max(c0 - 2 * w, 14.0)) + (0 if c0 > 20 else rvm(max(c0 - 2 * w, 14.0)) - rvm(14.0))
    Neff = rvm(c1 + 2 * w) - rvm(max(c0 - 2 * w, 2 * np.pi))
    print(f"RvM count in effective band ~ {Neff:.1f}; predicted Nyquist spacing ~ {(c1-c0)/Neff:.2f}")
    print(f"{'s':>5} {'n':>3} {'true':>10} {'gapR3':>10} {'gapAdapt':>10} {'R_adapt':>7} {'exact set (adaptive)':>30}")
    smax = 6.5
    for s in np.arange(smax, 1.19, -0.1):
        cs = np.arange(c0, c1 + 1e-9, s)
        if len(cs) < 4:
            continue
        f = L.build_forms(cs, w)
        A = f['archpole']; P = f['prime']
        true_m = float(np.linalg.eigvalsh(f['full'])[0])
        if true_m < 1e-13:
            print(f"{s:5.1f} {len(cs):3d} {true_m:+10.3e}  (past rank limit — stop)")
            break
        terms = prime_terms(cs)
        order = sorted(terms, key=lambda n: -np.linalg.norm(terms[n], 2))
        # R=3 chain
        Pres3 = sum(terms[n] for n in order[:3]); c_r3 = gap_bound(A - Pres3, P - Pres3, k=3)
        # adaptive: grow R until certified or R=12
        best = None
        for R in range(1, min(13, len(order) + 1)):
            Pres = sum(terms[n] for n in order[:R])
            b = gap_bound(A - Pres, P - Pres, k=3)
            if best is None or b > best[1]:
                best = (R, b, order[:R])
            if b > 0:
                break
        Rr, br, setr = best
        print(f"{s:5.1f} {len(cs):3d} {true_m:+10.3e} {c_r3:+10.3e} {br:+10.3e} {Rr:7d} {str(sorted(setr)):>30}",
              flush=True)
