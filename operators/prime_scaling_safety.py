#!/usr/bin/env python3
"""§856: the PRIME-SCALING SAFETY MARGIN t*. Find the largest t such that arch+pole − t·prime ⪰ 0 (all
generalized eigenvalues vs G ≥ 0). t=1 is the real Weil form. t* quantifies how close the ACTUAL primes
are to breaking Weil positivity: t*≫1 = primes comfortably inside positivity; t*→1 = on the boundary
(RH-marginal). Maps positivity-robustness in prime-configuration space — a new lens. Full Weil form,
valid (no §853 confound). Bisection on t."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam
NMAX = 3000; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

def mats(a, b, W, N, dr=0.004):
    fourW = 4*W
    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    CEN = list(np.linspace(a, b, N))
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    P = np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]
    cosu = np.cos(np.outer(ln, RGRID)); prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
    prime = 0.5*(prime+prime.T)
    def fc(c):
        y = W*(0.5j-c); y2 = W*(0.5j+c); return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in CEN]); pole = 2.0*np.real(np.outer(fp, fp))
    return pole+arch, prime, G

def tstar(AP, prime, G):
    # largest t with min gen-eig(AP - t*prime, G) >= 0
    def ok(t):
        return float(eigh(AP - t*prime, G, eigvals_only=True)[0]) >= 0
    if not ok(0.0): return 0.0
    lo, hi = 0.0, 1.0
    while ok(hi) and hi < 1e6: lo, hi = hi, hi*2
    for _ in range(50):
        mid = (lo+hi)/2
        if ok(mid): lo = mid
        else: hi = mid
    return lo

print("="*96)
print("§856 prime-scaling safety margin t* (largest t with arch+pole − t·prime ⪰ 0); t=1 is the real Weil form")
print("="*96)
print(f"{'window':>12} {'4W':>4} {'N':>3} {'N/Ncap':>7} {'margin(t=1)':>12} {'t* (safety)':>11}  interpretation")
for (a, b, fourW, frac) in [(14,44,2.0,0.7),(14,44,3.0,0.7),(14,44,4.0,0.7),(14,44,4.0,0.9),
                            (60,90,4.0,0.7),(120,150,4.0,0.7),(300,330,4.0,0.7)]:
    W = fourW/4; L = b-a; Ncap = L*2*W/np.pi; N = max(4, round(frac*Ncap))
    AP, prime, G = mats(a, b, W, N)
    m1 = float(eigh(AP - prime, G, eigvals_only=True)[0])
    ts = tstar(AP, prime, G)
    interp = "on boundary (RH-marginal)" if ts < 1.05 else ("comfortable" if ts < 3 else "very safe")
    print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {N:>3} {N/Ncap:>7.2f} {m1:>+12.4e} {ts:>11.4f}  {interp}")
print("="*96)
print("t* ≈ 1 → the actual primes sit right on the positivity boundary (RH-marginal, the hard corner).")
print("t* ≫ 1 → primes comfortably inside positivity (robust regime). Maps the safety margin across parameter space.")
print("="*96)
