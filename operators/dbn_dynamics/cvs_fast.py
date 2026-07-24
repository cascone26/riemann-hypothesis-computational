#!/usr/bin/env python3
"""
cvs_fast.py — reusable fast tooling for the CvS operator (LEDGER §170 tooling). SAFE speedups only:
 (#3) Precompute the t-independent arch+pole components ONCE per interval cL, and per-prime-power prime
      contributions once — so any prime-subset / coupling reassembly is cheap (exact reuse, no approximation).
 Instrumentation: reports the build (arch quadrature) vs solve (eigsy) time split so we KNOW the bottleneck
      before optimizing the solver. λ_min via mp.eigsy is kept as the TRUSTED path (robust to the near-null +
      near-degenerate tower, §167, and correct-signed unlike naive inverse iteration).
 (#2 drop-N, #4 cut-dps/T are LEFT OUT on purpose: silent-artifact risk — see LEDGER §169c, Groskin §6.6.)

API:
  comp = precompute_components(cL, N, T, dps)            # arch+pole+per-pp, timed  (the expensive part, once per cL)
  Q_even = assemble_even(comp, cP, N)                     # cheap: fold prime subset ≤cP into the even-sector matrix
  lam = lambda_min(Q_even)                                # trusted eigsy smallest eigenvalue
"""
import time
import mpmath as mp
import connes_cvs.operator as O

def _flint_prec(dps):
    if O.HAS_FLINT:
        from flint import ctx as _fctx
        _fctx.prec = int(dps * 3.5)

def precompute_components(cL, N, T, dps):
    """Compute (once) the arch+pole value/deriv per basis index and per-prime-power value/deriv.
    Prime contributions are keyed by prime-power n so any cP-subset is a cheap sum. Returns a dict."""
    mp.mp.dps = dps
    _flint_prec(dps)
    L = mp.log(cL)
    idx = list(range(-N, N + 1))
    arcpol = {}; arcpol_d = {}
    t0 = time.time()
    for n in idx:                                   # arch (expensive) + pole, t/prime-independent
        x = mp.mpf(n)
        O._hplus_cache_clear()
        arcpol[n]   = O.psi_pole(x, L) + O.psi_arch(x, L, T, dps)
        arcpol_d[n] = O.psi_pole_deriv(x, L) + O.psi_arch_deriv(x, L, T, dps)
        O._hplus_cache_clear()
    t_arch = time.time() - t0
    # per-prime-power contributions on this interval (cheap): keyed by pp n, using ALL pp up to a large cap.
    # caller picks cP; we compute for the full prime set ≤ int(cL) (superset of any cP≤cL of interest).
    pdata, _ = O.prime_powers_up_to(int(cL))
    pp = {}
    for (nn, logn, w) in pdata:
        val = {}; der = {}
        for n in idx:
            val[n] = O.psi_prime(mp.mpf(n), L, [(nn, logn, w)])
            der[n] = O.psi_prime_deriv(mp.mpf(n), L, [(nn, logn, w)])
        pp[nn] = (val, der)
    return {"cL": cL, "N": N, "L": L, "idx": idx,
            "arcpol": arcpol, "arcpol_d": arcpol_d, "pp": pp, "t_arch": t_arch}

def assemble_even(comp, cP, N):
    """Assemble the even-sector (N+1)x(N+1) matrix for prime set ≤ cP, reusing precomputed components."""
    idx = comp["idx"]; arcpol = comp["arcpol"]; arcpol_d = comp["arcpol_d"]; pp = comp["pp"]
    keep = [nn for nn in pp if nn <= cP]
    psi = {}; psid = {}
    for n in idx:
        pv = sum((pp[nn][0][n] for nn in keep), mp.mpf(0))
        pd = sum((pp[nn][1][n] for nn in keep), mp.mpf(0))
        psi[n]  = pv + arcpol[n]
        psid[n] = pd + arcpol_d[n]
    DIM = 2 * N + 1
    Q = mp.matrix(DIM, DIM)
    for i in range(DIM):
        m = i - N
        for j in range(DIM):
            nn = j - N
            Q[i, j] = psid[nn] if m == nn else (psi[m] - psi[nn]) / (m - nn)
    for i in range(DIM):
        for j in range(i + 1, DIM):
            a = (Q[i, j] + Q[j, i]) / 2; Q[i, j] = a; Q[j, i] = a
    # even-sector projection
    V = mp.matrix(DIM, N + 1); V[N, 0] = mp.mpf(1); r2 = 1 / mp.sqrt(2)
    for k in range(1, N + 1):
        V[N + k, k] = r2; V[N - k, k] = r2
    return V.T * Q * V

def lambda_min(Q_even):
    """Trusted smallest eigenvalue via mpmath symmetric eigensolver (robust to near-null + degenerate tower)."""
    ev = mp.eigsy(Q_even, eigvals_only=True)
    return min(ev[i] for i in range(Q_even.rows))

def lambda_min_gated(Q_even, tol=None):
    """Smallest eigenvalue WITH a convergence gate (LEDGER §170). Deep in the near-degenerate tower the
    matrix entries (adaptive mp.quad) are only accurate to the quadrature tolerance, so a reported λ_min far
    below that is unreliable — every §170 failure was eigsy returning a too-SHALLOW cluster member.

    Returns (lam, residual, trust): residual = ‖Q v − λ v‖ / ‖Q‖ for the smallest eigenPAIR; trust=False when
    residual > tol (default 10^(−dps/3), a conservative proxy for the quad accuracy floor). A False means the
    depth is NOT trustworthy — re-run at higher dps / tighter quadrature, do not report the number as-is."""
    n = Q_even.rows
    ev, V = mp.eigsy(Q_even)
    k = min(range(n), key=lambda i: ev[i])          # index of the algebraically smallest eigenvalue
    lam = ev[k]
    v = mp.matrix([V[r, k] for r in range(n)])
    Qv = Q_even * v
    res_vec = Qv - lam * v
    rnorm = mp.sqrt(sum(res_vec[r] ** 2 for r in range(n)))
    Qnorm = mp.sqrt(sum(Q_even[i, j] ** 2 for i in range(n) for j in range(n)))
    residual = rnorm / Qnorm if Qnorm != 0 else rnorm
    if tol is None:
        tol = mp.mpf(10) ** (-(mp.mp.dps // 3))
    return lam, residual, (residual <= tol)

# NOTE (§170): the residual gate catches eigsy NON-CONVERGENCE only. It does NOT catch the deeper failure —
# matrix entries inaccurate below the adaptive-quadrature tolerance — because eigsy converges fine on the
# (wrong) matrix, giving a small residual on a wrong λ. The DEFINITIVE trust check is CROSS-RESOLUTION
# AGREEMENT: recompute at (N+Δ, dps+Δ) and trust λ_min only where the two agree (this is what caught the
# §169c c=21 artifact: N=100 gave −54, N=120 gave −104 → disagreement → the N=100 value was garbage).
def trust_by_agreement(cL, cP, N, T, dps, dN=20, dprec=50):
    """Definitive deep-λ_min trust check: recompute at (N, N+dN)/(dps, dps+dprec); return both + agreement flag."""
    c1 = precompute_components(cL, N, T, dps);        lam1 = lambda_min(assemble_even(c1, cP, N))
    c2 = precompute_components(cL, N + dN, T, dps + dprec); lam2 = lambda_min(assemble_even(c2, cP, N + dN))
    l1 = float(mp.log10(abs(lam1))); l2 = float(mp.log10(abs(lam2)))
    return lam1, lam2, abs(l1 - l2) < 1.0   # agree if log10|λ| within 1 order

if __name__ == "__main__":
    import sys
    cL = float(sys.argv[1]) if len(sys.argv) > 1 else 13
    cP = float(sys.argv[2]) if len(sys.argv) > 2 else 13
    N  = int(sys.argv[3]) if len(sys.argv) > 3 else 30
    dps= int(sys.argv[4]) if len(sys.argv) > 4 else 60
    T  = int(sys.argv[5]) if len(sys.argv) > 5 else 150
    cLi = int(cL) if float(cL).is_integer() else cL
    t0 = time.time(); comp = precompute_components(cLi, N, T, int(dps)); t_build = time.time() - t0
    t1 = time.time(); Qe = assemble_even(comp, int(cP), N); lam = lambda_min(Qe); t_solve = time.time() - t1
    print(f"cL={cLi} cP={int(cP)} N={N} dps={dps} T={T}: log10|lam|={float(mp.log10(abs(lam))):.4f}")
    print(f"  TIMING: build(arch+pole+pp)={t_build:.1f}s (arch={comp['t_arch']:.1f}s)  "
          f"assemble+eigsy={t_solve:.1f}s  -> build is {100*t_build/(t_build+t_solve):.0f}% of total")
