#!/usr/bin/env python3
"""
cvs_interval.py — LEDGER §169: DECOUPLE interval-length from prime-content to dissect the CvS convergence
mechanism (Groskin finding (e): "interval extension, not prime identity, drives convergence", r=−0.96).

Groskin had ONE clean decoupled datum: c=13→14 improves γ1 by 5.8 OOM with NO new prime (14 not prime),
purely from L=log c growing. His c-sweep can't separate L from the prime set further (c locks both). OUR
operator exposes the pieces, so we build with interval L=log(c_L) but primes ≤ c_P INDEPENDENTLY, and trace
whole controlled curves:
  (A) FIX primes ≤ c_P0, EXTEND interval c_L = c_P0, c_P0+1, ... : does γ1 error keep dropping from L alone?
  (B) FIX interval c_L0, ADD primes c_P = ..., c_L0 : does adding primes at fixed interval help?
If (A) converges and (B) doesn't ⟹ convergence is GEOMETRIC (interval/scaling limit), not arithmetic —
reframing the open convergence-proof (⟺ RH by Hurwitz) as a geometric problem.

NOTE (honest): decoupled (c_L≠c_P) is a NON-standard chimera operator; the prime frequencies (1−log n/L) and
pole/arch all use L=log(c_L). Principled special case: c_L>c_P with c_P the last prime = exactly Groskin's
non-prime-c operators (e.g. c=14 = primes≤13 on interval log14). Interpret the trend, cautiously.

REQUIRES near-c-limited N (N≳100, dps≳200 for c up to ~20) — else γ1 error is N-limited (~1e-17 at N=60) and
the interval effect is masked. Built for Olympus. Usage: cvs_interval.py MODE cL cP N dps T
  MODE=one  : single (cL,cP) point — prints λ_min + γ1 error.
"""
import sys
import mpmath as mp
sys.path.insert(0, ".venv_cvs/lib/python3.14/site-packages")
try:
    import connes_cvs.operator as O
except ModuleNotFoundError:
    # Olympus fallback: package installed in a sibling venv path
    import importlib, glob
    for p in glob.glob("/home/cobo/**/connes_cvs", recursive=True):
        sys.path.insert(0, p.rsplit("/connes_cvs", 1)[0]); break
    import connes_cvs.operator as O

def build_decoupled(cL, cP, N, T, dps):
    """CvS Galerkin matrix with interval L=log(cL) but prime set ≤ cP (decoupled)."""
    mp.mp.dps = dps
    if O.HAS_FLINT:                              # MUST match build_galerkin: flint h_plus needs the prec set
        from flint import ctx as _fctx
        _fctx.prec = int(dps * 3.5)
    L = mp.log(cL)
    prime_data, _ = O.prime_powers_up_to(int(cP))
    idx = list(range(-N, N + 1))
    psi = {}; psid = {}
    for n in idx:
        x = mp.mpf(n)
        O._hplus_cache_clear()
        psi[n]  = O.psi_prime(x, L, prime_data) + O.psi_pole(x, L) + O.psi_arch(x, L, T, dps)
        psid[n] = O.psi_prime_deriv(x, L, prime_data) + O.psi_pole_deriv(x, L) + O.psi_arch_deriv(x, L, T, dps)
        O._hplus_cache_clear()
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
    return Q, L

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "one"
    cL  = float(sys.argv[2]) if len(sys.argv) > 2 else 14
    cP  = float(sys.argv[3]) if len(sys.argv) > 3 else 13
    N   = int(sys.argv[4]) if len(sys.argv) > 4 else 40
    dps = int(sys.argv[5]) if len(sys.argv) > 5 else 80
    T   = int(sys.argv[6]) if len(sys.argv) > 6 else 300
    cLi = int(cL) if float(cL).is_integer() else cL
    cPi = int(cP) if float(cP).is_integer() else cP
    Q, L = build_decoupled(cLi, cPi, N, T, dps)
    lam, vec = O.compute_ground_state(Q)
    zs = O.extract_zeros(vec, float(L), n_zeros=1, dps=dps)
    e1 = zs[0]['error']
    g1err = mp.nstr(e1, 6) if e1 is not None else "None(unresolved)"
    print(f"cL={cLi} (L={float(L):.5f}) cP={cPi} N={N} dps={dps} T={T}: "
          f"lambda_min={mp.nstr(lam,6)}  γ1_error={g1err}", flush=True)
