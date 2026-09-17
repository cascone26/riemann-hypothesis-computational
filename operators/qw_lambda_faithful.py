#!/usr/bin/env python3
"""cy372/373 — FAITHFUL build of the CM Weil quadratic form QW_λ in the V_n basis, directly from
ccm_2511.22755 §671-711 (NOT the abstract Jacobi PW_λ of rt6_prolate_matrix.py — this is the actual
Weil form whose simple-even property is CM's open 'missing step' §1385).

QW_λ(V_n,V_m) = W02 + W_p + W_R  on modes n,m ∈ {-N..N}.
Basis U_n on [0,L], L = 2 log λ.  Even kernel q(U_n,U_m)(y), y∈[0,L] (Lemma polarize0):
   n≠m:  [sin(2πmy/L) - sin(2πny/L)] / (π(n-m))
   n=m:  2(L-y) cos(2πny/L) / L
   q(.)(0) = 2·δ_{nm}.
W02 = 32 L sinh²(L/4) (L² - 16π² m n) / [(L²+16π² m²)(L²+16π² n²)]        (rank-2)
W_p = Σ_{1<k≤e^L} Λ(k) k^{-1/2} q(U_n,U_m)(log k)                         (von Mangoldt / primes)
W_R = (ω(0)/2)(γ + log(4π (e^L-1)/(e^L+1))) + ∫_0^L [e^{x/2}ω(x)-ω(0)]/(e^x-e^{-x}) dx
GOAL: is the smallest eigenvalue simple (gap>0) and its eigenvector EVEN (v_{-n}=v_n)?  = CM missing step 1.
"""
import mpmath as mp
mp.mp.dps = 30
from sympy import primerange, factorint
import numpy as np

def build_QW(N, lam, verbose=True):
    L = 2*mp.log(lam)
    pi = mp.pi
    idx = list(range(-N, N+1)); D = len(idx)
    # --- even kernel q(U_n,U_m)(y) on [0,L] ---
    def q(n, m, y):
        if n == m:
            return 2*(L-y)*mp.cos(2*pi*n*y/L)/L
        return (mp.sin(2*pi*m*y/L) - mp.sin(2*pi*n*y/L))/(pi*(n-m))
    def om0(n, m):   # q(.)(0)
        return mp.mpf(2) if n == m else mp.mpf(0)
    # --- W02 (rank-2) ---
    def W02(n, m):
        num = 32*L*mp.sinh(L/4)**2*(L**2 - 16*pi**2*m*n)
        den = (L**2 + 16*pi**2*m**2)*(L**2 + 16*pi**2*n**2)
        return num/den
    # --- prime powers k with Λ(k), 1<k≤e^L ---
    Kmax = int(mp.floor(mp.e**L))
    prime_powers = []  # (k, Lambda(k)=log p, log k)
    for p in primerange(2, Kmax+1):
        pk = p
        while pk <= Kmax:
            prime_powers.append((pk, mp.log(p), mp.log(pk)))
            pk *= p
    def W_p(n, m):
        s = mp.mpf(0)
        for k, Lam, logk in prime_powers:
            s += Lam*mp.power(k, mp.mpf(-1)/2)*q(n, m, logk)
        return s
    # --- W_R (archimedean), 1D numeric integral (removable sing. at 0) ---
    gamma_E = mp.euler
    Rconst_log = mp.log(4*pi*(mp.e**L - 1)/(mp.e**L + 1))
    def W_R(n, m):
        w0 = om0(n, m)
        def integrand(x):
            wx = q(n, m, x)
            den = mp.e**x - mp.e**(-x)
            return (mp.e**(x/2)*wx - w0)/den
        # integrand has a REMOVABLE singularity at 0 (num→0 like the denom); split at a tiny offset
        # far below any real ε_λ so the truncated ∫_0^eps ≈ eps·O(1) does NOT floor the spectrum.
        integ = mp.quad(integrand, [mp.mpf('1e-30'), L])
        return (w0/2)*(gamma_E + Rconst_log) + integ

    Q = mp.zeros(D, D)
    for a, n in enumerate(idx):
        for b, m in enumerate(idx):
            if b < a:
                Q[a, b] = Q[b, a]; continue
            Q[a, b] = W02(n, m) - W_R(n, m) - W_p(n, m)   # CM eq (bombtest): Ψ = W02 - W_R - Σ_p W_p
    if verbose:
        print(f"  built QW_λ: N={N} λ={float(lam):.4f} L={float(L):.4f} dim={D} #primepowers={len(prime_powers)}")
    return Q, idx, float(L)

def analyze(Q, idx, L):
    D = len(idx)
    A = np.array([[float(Q[i, j]) for j in range(D)] for i in range(D)])
    A = (A + A.T)/2
    ev, evec = np.linalg.eigh(A)
    eps1, eps2 = ev[0], ev[1]
    v = evec[:, 0]
    # parity: even means v_{-n} == v_n (index reversal symmetric)
    vrev = v[::-1]
    # fix sign
    if np.dot(v, vrev) < 0: vrev = -vrev
    even_part = (v + v[::-1])/2; odd_part = (v - v[::-1])/2
    parity_even_frac = np.linalg.norm(even_part)/np.linalg.norm(v)
    print(f"  ε1={eps1:+.6e}  ε2={eps2:+.6e}  gap={eps2-eps1:.3e}  (simple? {eps2-eps1>1e-9})")
    print(f"  ground-state even-fraction = {parity_even_frac:.6f}  (1.0 = perfectly even)  → EVEN? {parity_even_frac>0.999}")
    return eps1, eps2, parity_even_frac

if __name__ == "__main__":
    print("=== FAITHFUL QW_λ (CM §671-711): simple-even check across λ,N ===")
    for lam, N in [(mp.e**2, 20), (mp.e**3, 24)]:
        Q, idx, L = build_QW(N, lam)
        analyze(Q, idx, L)
        print()
