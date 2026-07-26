#!/usr/bin/env python3
"""ccm_discrepancy.py -- §257: CCM 2310.18423 discrepancy coefficients d_n, gate-verified.

CCM cyclic pair (D=scaling op, xi=Gaussian). Spectral measure of the SEMILOCAL prolate program:
   dmu_S(s) ~ |L_inf(1/2+is) * prod_{p in S} (1-p^{-1/2-is})^{-1}|^2 ds
where |L_inf|^2 ~ |Gamma(1/4+is/2)|^2 (archimedean Euler factor, CCM line 216).
The arithmetic content = the DISCREPANCY d_n = beta_n(S) - beta_n(arch), beta_n(arch)=(2n-1)(2n)/4
(the monic Jacobi recurrence coeff; CCM a_n=1/2 sqrt((2n+1)(2n+2)), Thm dzero: d_n=0 iff S empty).

Method: modified Chebyshev / Wheeler with the ARCHIMEDEAN orthogonal polynomials as reference
(numerically stable to high n, unlike direct a_n^2 Stieltjes which loses digits by n~5).
GATE (must hold or discard): odd modified moments = 0 exactly (even measure) -- hard-zeroed, NOT
integrated (the 2*int_0^inf even-trick is INVALID for odd integrands; that bug gave a false arch d_n).
Arch d_n then come out ~1e-36..1e-59 (machine zero) = CCM Thm dzero reproduced + method validated.

RESULT (dps=60): low-n d_n match the independent direct-Stieltjes run (§257). d_n do NOT decay --
they oscillate in sign and stay O(1) relative to b_arch(n) (|d_n/b_arch(n)| ~ 0.5-0.7): the semilocal
deformation is an O(1) log-p-periodic MODULATION of the Jacobi coefficients, not a c_0(N)-vanishing
perturbation. That modulation structure is a concrete input to CCM's sought semilocal prolate operator.
"""
import sys, mpmath as mp
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 60
PI = mp.pi

def b_arch(k):
    return mp.mpf(1) if k == 0 else (2*k-1)*(2*k)/mp.mpf(4)

def measure(s, S):
    g = mp.gamma(mp.mpf(1)/4 + 1j*s/2); w = mp.fabs(g)**2
    for p in S:
        w *= 1/(1 + mp.mpf(1)/p - 2*p**(-mp.mpf(1)/2)*mp.cos(s*mp.log(p)))
    return w

def piarch(s, K):
    v = [mp.mpf(1)]
    if K >= 1: v.append(s)
    for k in range(1, K): v.append(s*v[k] - b_arch(k)*v[k-1])
    return v

def modmoments(S, K):
    per = [2*PI/mp.log(p) for p in S] if S else [10]
    nodes = [0]; x = mp.mpf(0); step = min(per)/5 if S else mp.mpf(2)
    while x < 90: x += step; nodes.append(x)
    nu = []
    for k in range(2*K):
        if k % 2 == 1: nu.append(mp.mpf(0))                 # GATE: odd = 0 by symmetry
        else: nu.append(2*mp.quad(lambda s: piarch(s, k)[k]*measure(s, S), nodes))
    return nu

def wheeler(nu, K):
    sig = {}
    for l in range(-1, 2*K): sig[(-1, l)] = mp.mpf(0)
    for l in range(2*K): sig[(0, l)] = nu[l]
    alpha = [nu[1]/nu[0]]; beta = [nu[0]]
    for k in range(1, K):
        for l in range(k, 2*K-k):
            sig[(k, l)] = (sig[(k-1, l+1)] - alpha[k-1]*sig[(k-1, l)]
                           - beta[k-1]*sig.get((k-2, l), mp.mpf(0)) + b_arch(l)*sig[(k-1, l-1)])
        alpha.append(sig[(k, k+1)]/sig[(k, k)] - sig[(k-1, k)]/sig[(k-1, k-1)])
        beta.append(sig[(k, k)]/sig[(k-1, k-1)])
    return alpha, beta

def d_coeffs(S, K=10):
    al, be = wheeler(modmoments(S, K), K)
    return [be[n] - b_arch(n) for n in range(1, K)]

if __name__ == "__main__":
    K = 10
    print(f"CCM discrepancy d_n (dps={mp.mp.dps}), n=1..{K-1}. arch row = GATE (must be ~0):")
    for S in [[], [2], [3], [2, 3], [2, 3, 5, 7]]:
        dn = d_coeffs(S, K)
        tag = "arch(GATE)" if not S else f"S={S}"
        print(f"  {tag:14s}: {[mp.nstr(x, 5) for x in dn]}")

# ============================================================
# CLOSED FORM (§261a-e, verified machine-precision):
#   nu_k = Sum_n a_n * pi_k^MP(-i d/dw)[ 2*sqrt(2)*pi^(3/2) * cosh(w)^(-1/2) ] |_{w=log n}
#   - Phi_0(w) = 2*sqrt(2)*pi^(3/2) / sqrt(cosh(w))   [FT of |Gamma(1/4+is/2)|^2, the arch measure; POSITIVE-DEFINITE]
#   - Phi_k(w) = pi_k^MP(-i d/dw) Phi_0(w)            [MP-polynomial as differential operator on sech^{1/2}]
#   - frequencies {log n} (integer lattice), Euler weights a_n = |zeta-partial(1/2+is)|^2 coefficients
#   Arch cyclic pair = Meixner-Pollaczek (lambda=1/4, phi=pi/2) EXACTLY (b_arch(n)=(2n-1)(2n)/4).
#   CCM<->Weil gap (RH locus): CCM = |zeta-Euler|^2 level (positive); Weil = von Mangoldt = zeta'/zeta level.
# ============================================================
