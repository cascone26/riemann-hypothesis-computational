#!/usr/bin/env python3
"""§969 GUARD: is κ = pole ORDER, or pole LOCATION-COUNT? test ζ, ζ², ζ³ (poles of order 1,2,3).  (HANDLE)

§958 verified κ(arch−prime)=pole order for orders {0,1} (entire χ vs ζ). ζ^k has a pole of ORDER k at s=1.
But −(ζ^k)'/ζ^k = k·ΣΛ (prime coeffs ×k) and the completed ξ^k has Γ-factor^k (arch density ×k), so
arch−prime for ζ^k = k·(arch−prime for ζ) — a POSITIVE SCALING, inertia-invariant ⇒ predict κ=1 for ALL k,
NOT k. If confirmed, the honest invariant is 'κ = number of pole LOCATIONS at s=1' (0 or 1 for everything in
the Selberg class, since a pole only ever sits at s=1), and §958's 'κ=pole order' held only because the tested
orders were {0,1}. Verify-before-ship the scaling claim rather than assert it. HANDLE (own code)."""
import numpy as np
from scipy.special import digamma

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam
NMAX = 4000; LAM = vonmangoldt(NMAX)
NS = np.arange(2, NMAX+1); LOGN = np.log(NS); LAMv = LAM[2:]

def form(power, cen_lo, cen_hi, fourW, frac=0.8, dr=0.005):
    """arch−prime for ζ^power: arch density = power·Ω, prime coeffs ×power."""
    W = fourW/4; L = cen_hi-cen_lo; N = max(4, round(frac*L*2*W/np.pi))
    RG = np.arange(-(cen_hi+40), (cen_hi+40), dr)
    OMEGA = power*(np.real(digamma(0.25 + 0.5j*RG)) - np.log(np.pi))
    CEN = np.linspace(cen_lo, cen_hi, N)
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    P = np.array([fej(RG-c)+fej(RG+c) for c in CEN]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; wts = power*LAMv[mask]/np.sqrt(NS[mask])
    cosu = np.cos(np.outer(ln, RG)); prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*wts[None, :]).sum(axis=1)
    return arch - 0.5*(prime+prime.T), G

def kappa(M, G):
    gw, gv = np.linalg.eigh(G); keep = gw > gw.max()*1e-12
    B = gv[:, keep]/np.sqrt(gw[keep]); ev = np.linalg.eigvalsh(B.T @ M @ B)
    s = max(1.0, np.abs(ev).max()); return int((ev < -1e-9*s).sum())

print("="*80)
print("§969 κ(arch−prime) for ζ^k. pole ORDER=k. 'order' predicts κ=k; 'location' predicts κ=1.")
print("="*80)
print(f"{'ζ^k':>6} {'pole order':>11} | {'κ @4W=5':>8} {'κ @4W=6':>8} {'κ @4W=7':>8}")
print("-"*80)
for k in [1, 2, 3]:
    ks = []
    for fw in [5.0, 6.0, 7.0]:
        M, G = form(k, 14, 44, fw); ks.append(kappa(M, G))
    print(f"{'ζ^'+str(k):>6} {k:>11} | {ks[0]:>8} {ks[1]:>8} {ks[2]:>8}")
print("-"*80)
print("READ: κ=1 for ALL k confirms arch−prime for ζ^k = k·(ζ form) (positive scaling, inertia-invariant) —")
print("so the invariant is κ = #pole LOCATIONS at s=1 (∈{0,1}), NOT pole order. §958's 'κ=pole order' held")
print("only for orders {0,1}. Honest refinement of §958b, guarding against overclaiming order-2+.")
print("="*80)
