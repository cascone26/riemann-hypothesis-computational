#!/usr/bin/env python3
"""§178e diagnostic: WHY does the bridge overlap 'starve' at deep c / high N?
Hypothesis: NOT roundoff — MODE-CROSSING. Prolate theory: A_λ has ~4c small eigenvalues (c=6→~24). At high N the
globally-smallest even eigenvalue may no longer be the ε(λ) mode that k_λ (an h_0/h_4 combo) approximates; cvs_ground
picks argmin(eigenvalue) → grabs the wrong mode → overlap collapses to ~0 ('starved').
This tool dumps the smallest few even eigenvalues AND each one's overlap+dist with k_λ, so we can SEE which eigenvalue
k_λ actually matches and where it sits in the ordering. If a k_λ-matching mode (overlap→1) exists but isn't the min,
the fix is to SELECT by k_λ-overlap (or by continuity), not to raise dps. Reuses the exact cvs_ground matrix build."""
import sys, time
import numpy as np
import mpmath as mp
sys.path.insert(0, ".venv_cvs/lib/python3.14/site-packages")
import connes_cvs.operator as O
from bridge_overlap import k_lambda_on, f_psi, overlap

def even_spectrum(c, N, T, dps):
    """Return (L, eigs_list, vecs_matrix) for the even-sector Weil form Qe — full spectrum, ascending eigsy order."""
    mp.mp.dps = dps
    if O.HAS_FLINT:
        from flint import ctx as _fctx; _fctx.prec = int(dps * 3.5)
    L = mp.log(c); pdata, _ = O.prime_powers_up_to(int(c))
    idx = list(range(-N, N + 1)); psiv = {}; psid = {}
    for m in idx:
        x = mp.mpf(m); O._hplus_cache_clear()
        psiv[m] = O.psi_prime(x, L, pdata) + O.psi_pole(x, L) + O.psi_arch(x, L, T, dps)
        psid[m] = O.psi_prime_deriv(x, L, pdata) + O.psi_pole_deriv(x, L) + O.psi_arch_deriv(x, L, T, dps)
        O._hplus_cache_clear()
    DIM = 2 * N + 1; Q = mp.matrix(DIM, DIM)
    for i in range(DIM):
        m = i - N
        for j in range(DIM):
            nn = j - N
            Q[i, j] = psid[nn] if m == nn else (psiv[m] - psiv[nn]) / (m - nn)
    for i in range(DIM):
        for j in range(i + 1, DIM):
            a = (Q[i, j] + Q[j, i]) / 2; Q[i, j] = a; Q[j, i] = a
    Vv = mp.matrix(DIM, N + 1); Vv[N, 0] = mp.mpf(1); r2 = 1 / mp.sqrt(2)
    for k in range(1, N + 1): Vv[N + k, k] = r2; Vv[N - k, k] = r2
    Qe = Vv.T * Q * Vv
    eigs, vecs = mp.eigsy(Qe)
    return float(L), eigs, vecs

if __name__ == "__main__":
    c   = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    N   = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    dps = int(sys.argv[3]) if len(sys.argv) > 3 else 90
    T   = int(sys.argv[4]) if len(sys.argv) > 4 else 400
    K   = int(sys.argv[5]) if len(sys.argv) > 5 else 8      # how many smallest even modes to inspect
    M   = 12000
    lam = float(np.sqrt(c))
    t0 = time.time()
    L, eigs, vecs = even_spectrum(c, N, T, dps)
    order = sorted(range(N + 1), key=lambda i: float(eigs[i]))   # ascending eigenvalue
    x = np.linspace(0.0, L, M + 1)
    kA = k_lambda_on(x - L / 2.0, lam)
    print(f"§178e diag  c={c} N={N} dps={dps} T={T}  L=log c={L:.4f}  λ={lam:.4f}  (~{4*c} small eigs expected)  [{time.time()-t0:.0f}s]")
    print(f"{'rank':>4} {'log10 eig':>10} {'overlap(k_λ)':>13} {'dist(k_λ)':>11}   {'← k_λ mode?':>12}")
    best = (-1, -1.0, None)
    for r in range(min(K, N + 1)):
        i = order[r]
        ev = [float(vecs[row, i]) for row in range(N + 1)]
        fp = f_psi(x, ev, L)
        ov, d = overlap(fp, kA, x)
        loge = float(mp.log10(abs(eigs[i])))
        flag = "  ***" if ov > 0.99 else ""
        if ov > best[1]: best = (r, ov, (loge, d))
        print(f"{r:>4} {loge:>10.2f} {ov:>13.6f} {d:>11.3e}{flag}")
    print(f"\nBEST k_λ-match: rank {best[0]} in the ascending spectrum, overlap {best[1]:.6f}, "
          f"log10 eig {best[2][0]:.2f}, dist {best[2][1]:.3e}")
    if best[0] == 0:
        print("→ k_λ IS the global minimum here: NOT mode-crossing at this (c,N). argmin selector is correct.")
    else:
        print(f"→ MODE-CROSSING CONFIRMED: k_λ matches rank {best[0]}, NOT the global min. argmin grabs the wrong mode.")
        print("   FIX for §178e: select the eigenvector by max overlap with k_λ (or track by continuity), not argmin.")
