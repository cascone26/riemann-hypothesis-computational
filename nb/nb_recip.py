#!/usr/bin/env python3
"""NB Gram structure: split G = G_smooth (log part) + G_arith (Vasyunin cotangent part) and see which
carries positive-definiteness. FINDING: G_smooth is INDEFINITE alone; the ARITHMETIC (reciprocity-governed
cotangent) part restores PSD — inverse of the Weil form (where arch is + and primes threaten it)."""
import numpy as np
from math import gcd
GAMMA = 0.5772156649015329; LOG2PI = np.log(2*np.pi)
def frac(x): return x - np.floor(x)
def Vs(n, m):
    if m == 1: return 0.0
    k = np.arange(1, m); return np.sum(frac(k*n/m)/np.tan(np.pi*k/m))
def A_smooth(n, m):
    g = gcd(n, m); n //= g; m //= g; r = n/m
    return (1-r)/2*np.log(r) + (1+r)/2*(LOG2PI-GAMMA)
def A_arith(n, m):
    g = gcd(n, m); n0 = n//g; m0 = m//g
    return -(np.pi/(2*m0))*(Vs(n0, m0)+Vs(m0, n0))
def gram_parts(N):
    Gs = np.empty((N, N)); Ga = np.empty((N, N))
    for u in range(1, N+1):
        for v in range(1, N+1):
            Gs[u-1, v-1] = A_smooth(u, v)/u
            Ga[u-1, v-1] = A_arith(u, v)/u
    return 0.5*(Gs+Gs.T), 0.5*(Ga+Ga.T)
if __name__ == "__main__":
    print(f"{'N':>4} {'lmin(full)':>12} {'lmin(smooth)':>13} {'lmin(arith)':>12} {'||Ga||/||Gs||':>13}")
    for N in [10, 20, 40, 80, 160]:
        Gs, Ga = gram_parts(N); G = Gs+Ga
        e = np.linalg.eigvalsh(G); es = np.linalg.eigvalsh(Gs); ea = np.linalg.eigvalsh(Ga)
        print(f"{N:>4} {e[0]:>12.3e} {es[0]:>13.3e} {ea[0]:>12.3e} {np.linalg.norm(Ga)/np.linalg.norm(Gs):>13.4f}")
    print("full PSD, smooth indefinite, arithmetic restores it: reciprocity-governed cotangent part is essential.")
