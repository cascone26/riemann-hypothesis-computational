#!/usr/bin/env python3
"""FUNCTION-FIELD ROSETTA, instrument 1: the exact Weil form for elliptic curves over F_p.
RH here is PROVEN (Hasse): zeta_E(T) = (1-aT+pT^2)/((1-T)(1-pT)), zeros T=p^{-1/2}e^{±i theta},
a = p+1-#E(F_p), |a| <= 2 sqrt(p). Frobenius eigenvalues alpha, bar-alpha = sqrt(p) e^{±i theta}.

Explicit formula (exact, finite): for test coefficients c_1..c_K on degrees n,
   #E(F_{p^n}) = p^n + 1 - alpha^n - bar-alpha^n
so    alpha^n + bar-alpha^n = 2 p^{n/2} cos(n theta) = p^n + 1 - #E(F_{p^n})  := a_n  (integers!).
The Weil positivity analog: for the completed zeta, the quadratic form
   Q(c) = sum_{n,m} c_n c_m [zero side]  with zero side = sum over the two zeros of ...
Normalized zero side (unitarized, zeros e^{±i theta} on the unit circle after T -> T/sqrt(p)):
   Q_zero(c) = sum_{zeros z=e^{±i theta}} |C(z)|^2,  C(z) = sum_n c_n z^n     -- PSD trivially, but
   the ARITHMETIC identity is: |C(e^{i th})|^2 + |C(e^{-i th})|^2 expressed via the a_n (point counts):
   sum_{n,m} c_n c_m (e^{i(n-m)th}+e^{-i(n-m)th}+e^{i(n+m)th}+e^{-i(n+m)th})... For REAL c and the pair,
   Q_zero(c) = 2 sum_{n,m} c_n c_m [cos((n-m)th) + cos((n+m)th)]  and  2cos(k th) = a_k / p^{k/2}   (a_0=2).
So the EXACT function-field Weil form is
   W[n,m] = ahat_{|n-m|} + ahat_{n+m},   ahat_k = a_k / p^{k/2}  (ahat_0 = 2),
and RH(=Hasse, proven) <=> W is PSD for ALL curves and all K  <=> |theta| real.
This mirrors the number-field 'full' form (cos(r log n) structure with r=zeros <-> here theta).
QUESTIONS the Rosetta asks: (1) verify W PSD for many curves (sanity); (2) MARGIN structure: min-eig
of W across the family a in [-2 sqrt p, 2 sqrt p] -- where is the razor (a=+-2 sqrt p supersingular
extremal = theta=0,pi: DEGENERATE zeros!); (3) what does the PROVEN mechanism (Hasse/Castelnuovo)
look like at instrument level: margin vs distance-to-extremal.
"""
import numpy as np

def frob_angle(a, p):
    x = a/(2*np.sqrt(p))
    assert abs(x) <= 1+1e-12, "violates Hasse -- impossible for a real curve"
    return np.arccos(np.clip(x, -1, 1))

def weil_form(a, p, K):
    th = frob_angle(a, p)
    ah = lambda k: 2*np.cos(k*th)            # = a_k / p^{k/2} exactly
    W = np.empty((K, K))
    for n in range(1, K+1):
        for m in range(1, K+1):
            W[n-1, m-1] = ah(abs(n-m)) + ah(n+m)
    return W

if __name__ == "__main__":
    import sys
    p = int(sys.argv[1]) if len(sys.argv) > 1 else 101
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    amax = int(np.floor(2*np.sqrt(p)))
    print(f"FF Weil form, p={p}, K={K}, Hasse interval a in [-{amax},{amax}] (all values realized by curves)")
    print(f"{'a':>5} {'theta/pi':>9} {'min-eig(W)':>12} {'PSD?':>5}")
    rows = []
    for a in range(-amax, amax+1):
        W = weil_form(a, p, K)
        ev = np.linalg.eigvalsh(W)
        rows.append((a, ev[0]))
        if abs(a) in (amax, amax-1, amax//2, 1, 0):
            print(f"{a:>5} {frob_angle(a,p)/np.pi:>9.4f} {ev[0]:>12.4e} {str(ev[0] > -1e-9):>5}")
    r