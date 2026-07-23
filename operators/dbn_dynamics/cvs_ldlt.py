#!/usr/bin/env python3
"""§179: is the §176 prime cascade the PIVOT structure of the CvS Weil matrix's LDL^T?
The cascade decomposes ONE near-null eigenvalue over primes (Hellmann-Feynman, an UPPER bound — variational).
LDL^T Q=LDL^T is manifestly PD iff every pivot d_i>0 — pivots are 'what survives after cancelling everything
before,' structurally the same 'telescoping' operation. If the pivot descent mirrors the cascade (smallest
pivot ≈ λmin, reached by a telescoping that matches the per-prime h(n) drops), the cascade becomes the
anatomy of a CONSTRUCTIVE positivity certificate (lower bound), not just an eigenvector story.
Runs on the SAME CvS even-sector matrix as cvs_gate/bridge (connes_cvs.operator), at cascade-adequate dps
(prior A1 tool was dps=60 for c=14 = noise-starved). Reports pivots vs known λmin and the Fuchs −4πc depth."""
import sys, time, math
import mpmath as mp
sys.path.insert(0, ".venv_cvs/lib/python3.14/site-packages")
import connes_cvs.operator as O

def build_cvs_even(c, N, T, dps):
    """Even-sector CvS Weil matrix Qe (same build as bridge_overlap.cvs_ground / cvs_gate)."""
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
    return Vv.T * Q * Vv, float(L)

def ldlt_pivots(Q, D):
    d = [mp.mpf(0)] * D; Lf = [[mp.mpf(0)] * D for _ in range(D)]
    for i in range(D):
        s = Q[i, i]
        for k in range(i): s -= Lf[i][k] * Lf[i][k] * d[k]
        d[i] = s
        for j in range(i + 1, D):
            t = Q[j, i]
            for k in range(i): t -= Lf[j][k] * Lf[i][k] * d[k]
            Lf[j][i] = t / d[i] if d[i] != 0 else mp.mpf(0)
    return d

if __name__ == "__main__":
    c   = int(sys.argv[1]) if len(sys.argv) > 1 else 13
    N   = int(sys.argv[2]) if len(sys.argv) > 2 else 80
    dps = int(sys.argv[3]) if len(sys.argv) > 3 else 90
    T   = int(sys.argv[4]) if len(sys.argv) > 4 else 400
    t0 = time.time()
    D = N + 1
    Qe, L = build_cvs_even(c, N, T, dps)
    d = ldlt_pivots(Qe, D)
    dv = [mp.log10(abs(x)) if x != 0 else mp.mpf('-999') for x in d]
    dvf = [float(x) for x in dv]
    sgn = [1 if x > 0 else (-1 if x < 0 else 0) for x in d]
    nneg = sum(1 for s in sgn if s < 0)
    fuchs = 6.37 - 4 * math.pi * c / math.log(10) + 4.5 * math.log(c) / math.log(10)  # predicted log10 λmin
    order = sorted(range(D), key=lambda i: dvf[i])
    print(f"§179 CvS LDL^T pivots  c={c} N={N} dps={dps} T={T}  D={D}  L={L:.4f}  [{time.time()-t0:.0f}s]")
    print(f"  #pivots<=0 : {nneg}   (must be 0 for PD / a positivity certificate)")
    print(f"  min log10|pivot| = {min(dvf):.3f}   Fuchs-predicted log10 λmin = {fuchs:.3f}   (match ⇒ smallest pivot IS the near-null)")
    print(f"  Σ log10|pivot| (=log10 det) = {sum(dvf):.3f}")
    print(f"  smallest 12 log10|pivot|: {[f'{dvf[i]:.2f}' for i in order[:12]]}")
    print(f"  #|pivot|<1e-12 of full: {sum(1 for x in dvf if x < -12)}   (near-kernel size — cascade is ONE mode, a big near-kernel would mean NO single-pivot cascade)")
    print(f"\n  DIAGONAL-ORDER pivot descent (does it telescope like the §176 cascade?):")
    print(f"  {'i':>4} {'log10|d_i|':>11} {'sign':>5}")
    for i in list(range(min(8, D))) + ['...'] + list(range(D - 6, D)):
        if i == '...': print("   ...")
        else: print(f"  {i:>4} {dvf[i]:>11.3f} {sgn[i]:>5}")
    print("\nREAD: (1) all pivots>0 = independent positivity certificate at this c. (2) if min pivot ≈ Fuchs λmin AND")
    print("a SINGLE pivot carries it (small near-kernel), the near-null is one pivot ⇒ next: does that pivot's Schur")
    print("formation telescope per-prime like §176. (3) large near-kernel ⇒ cascade≠pivots, pivot to zero-density (§123).")
