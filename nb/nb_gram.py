#!/usr/bin/env python3
"""Nyman-Beurling / Baez-Duarte Gram matrix (q=1), Vasyunin closed form, VALIDATED against the direct
|zeta|^2 integral. Formulas locked from Bettin-Conrey (arXiv:1111.0931) + Ehm 2024 (arXiv:2405.06349).

Space: functions phi_n(s) = n^{-s} zeta(s)/s on Re s=1/2, inner product <f,g>=(1/2pi) int f g-bar dt.
Direct:  G[u,v] = (uv)^{-1/2} (1/2pi) int (u/v)^{it} |zeta(1/2+it)|^2 / (1/4+t^2) dt.
Vasyunin (coprime n,m): A(n/m) = (1-r)/2 log r + (1+r)/2 (log2pi-gamma) - (pi/2m)(V(n,m)+V(m,n)),
   r=n/m,  V(n,m) = sum_{k=1}^{m-1} frac(kn/m) cot(pi k/m).
Reciprocity check: G[u,v] symmetric <=> S_1(1/r)=r S_1(r)+K(1-r)+(1/2)(1+r)log r (baked in).
GOAL of this file: nail the exact relation between direct G[u,v] and Vasyunin A(n/m) (indexing/normalization),
so the fast closed-form Gram is trustworthy. Then nb_spectrum.py studies conditioning + the reciprocity control.
"""
import numpy as np
import mpmath as mp
mp.mp.dps = 25

GAMMA = float(mp.euler); LOG2PI = float(mp.log(2*mp.pi))
K1 = (LOG2PI - GAMMA + 1)/2

def frac(x): return x - np.floor(x)

def vasyunin_V(n, m):
    if m == 1: return 0.0
    k = np.arange(1, m)
    return np.sum(frac(k*n/m) * (1.0/np.tan(np.pi*k/m)))

def vasyunin_A(n, m):
    """Autocorrelation for coprime n,m (r=n/m)."""
    from math import gcd
    g = gcd(n, m); n//=g; m//=g
    r = n/m
    return (1-r)/2*np.log(r) + (1+r)/2*(LOG2PI-GAMMA) - (np.pi/(2*m))*(vasyunin_V(n,m)+vasyunin_V(m,n))

def direct_G(u, v, T=400.0):
    """(uv)^{-1/2}(1/2pi) int (u/v)^{it} |zeta(1/2+it)|^2/(1/4+t^2) dt, real part (integrand is even+Hermitian)."""
    lr = mp.log(mp.mpf(u)/v)
    def integrand(t):
        z = mp.zeta(mp.mpf('0.5') + 1j*t)
        return (mp.e**(1j*lr*t) * abs(z)**2 / (mp.mpf('0.25')+t*t)).real
    I = mp.quad(integrand, [-T, -50, -10, 0, 10, 50, T])
    return float((mp.mpf(u)*v)**mp.mpf('-0.5') * I/(2*mp.pi))

def G(u, v):
    """VALIDATED closed-form Nyman-Beurling-Baez-Duarte Gram entry (q=1): G[u,v] = A(u,v)/u.
    Symmetric via reciprocity: A(u,v)/u = A(v,u)/v. No zeta integrals."""
    return vasyunin_A(u, v)/u

def build_gram(N):
    M = np.empty((N, N))
    for u in range(1, N+1):
        for v in range(1, N+1):
            M[u-1, v-1] = G(u, v)
    return 0.5*(M+M.T)   # enforce exact symmetry (reciprocity holds analytically; guards float)

if __name__ == "__main__":
    import sys
    if "val" in sys.argv:
        print("VALIDATION: direct |zeta|^2 integral vs G[u,v]=A(u,v)/u")
        print(f"{'(u,v)':>10} {'direct':>12} {'A(u,v)/u':>12} {'rel.err':>10}")
        for (u,v) in [(1,1),(2,1),(2,3),(3,2),(3,4),(2,5),(5,6),(4,6),(6,4)]:
            d = direct_G(u,v); g = G(u,v)
            print(f"{str((u,v)):>10} {d:>12.6f} {g:>12.6f} {abs(d-g)/abs(g):>10.2e}")
    else:
        print("Nyman-Beurling-Baez-Duarte Gram spectrum vs N (validated closed form, no zeta integrals):")
        print(f"{'N':>4} {'lambda_min':>13} {'lambda_max':>13} {'cond':>12} {'det^(1/N)':>12}")
        for N in [5, 10, 20, 40, 80, 160]:
            M = build_gram(N)
            ev = np.linalg.eigvalsh(M)
            cond = ev[-1]/ev[0] if ev[0] > 0 else np.inf
            sign, logdet = np.linalg.slogdet(M)
            print(f"{N:>4} {ev[0]:>13.4e} {ev[-1]:>13.4e} {cond:>12.3e} {np.exp(logdet/N):>12.4e}", flush=True)
        print("\n(PSD [lambda_min>0] required — Gram of independent functions. Conditioning growth = the")
        print(" numerical wall the reciprocity structure might tame; det^(1/N)->0 rate ~ d_N behavior.)")
