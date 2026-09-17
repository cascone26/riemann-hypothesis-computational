#!/usr/bin/env python3
"""§964 THE κ≤1 PERSISTENCE TEST — does the ONE negative square stay pinned as bandwidth grows?

The whole burden of RH (in the §947-958 reframing) is κ(arch−prime)≤1. Every static probe at fixed
bandwidth keeps reducing to the §891 alignment razor. But κ≤1 is fundamentally a LIMIT statement: as the
band-limit W→∞ (more prime-power terms enter, the compression matrix grows), does #neg(arch−prime) STAY
EXACTLY 1, or does it climb? If it stays pinned at 1 across a wide, growing sweep, that is the strongest
DIRECT numerical evidence for the upper bound (and the scaling of the most-negative eigenvalue vs the
smallest-positive 'gap' exposes HOW the single negative direction stays isolated). If it climbs, the
'one negative square' is a low-bandwidth artifact and the framing needs the pole handled differently.

For ζ: build the Weil form arch−prime (no pole) on a window, sweep 4W and let N scale with it. Track:
  κ = #neg(arch−prime);  λ_min (the negative eig);  λ_+1 (smallest POSITIVE eig, the 'gap' above 0);
  ratio |λ_min|/λ_+1 (how deep the single negative sits vs the positive floor).
HANDLE (own code) — report whatever the numbers do, honestly, either way."""
import numpy as np
from scipy.special import digamma

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam
NMAX = 40000; LAM = vonmangoldt(NMAX)
NS = np.arange(2, NMAX+1); LOGN = np.log(NS); LAMv = LAM[2:]

def weil_form(cen_lo, cen_hi, fourW, frac, dr=0.004):
    """ζ arch−prime form (no pole). N scales with bandwidth*window."""
    W = fourW/4; L = cen_hi-cen_lo; Ncap = L*2*W/np.pi; N = max(4, round(frac*Ncap))
    RG = np.arange(-(cen_hi+120), (cen_hi+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RG)) - np.log(np.pi)         # ζ archimedean (q=1)
    CEN = np.linspace(cen_lo, cen_hi, N)
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    P = np.array([fej(RG-c)+fej(RG+c) for c in CEN]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; wts = LAMv[mask]/np.sqrt(NS[mask]); m = int(mask.sum())
    cosu = np.cos(np.outer(ln, RG)); prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*wts[None, :]).sum(axis=1)
    prime = 0.5*(prime+prime.T)
    return arch - prime, G, N, m

def inertia(M, G):
    gw, gv = np.linalg.eigh(G); keep = gw > gw.max()*1e-12
    B = gv[:, keep]/np.sqrt(gw[keep]); ev = np.sort(np.linalg.eigvalsh(B.T @ M @ B))
    s = max(1.0, np.abs(ev).max()); tol = 1e-9*s
    kneg = int(np.sum(ev < -tol)); kpos = int(np.sum(ev > tol)); kzero = len(ev)-kneg-kpos
    lam_min = ev[0]
    posvals = ev[ev > tol]; lam_p1 = posvals.min() if posvals.size else np.nan
    return kneg, kpos, kzero, lam_min, lam_p1, len(ev)

print("="*104)
print("§964 κ(arch−prime) PERSISTENCE vs bandwidth — does the ONE negative square stay pinned as W grows?")
print("window [14,44], N scales as 0.8·(b−a)·2W/π. ζ (q=1), no pole term.")
print("="*104)
print(f"{'4W':>5} {'#primes':>8} {'N':>4} {'dim':>4} | {'κ=#neg':>7} {'#pos':>5} {'#zero':>6} | "
      f"{'λ_min':>11} {'λ_+1 (gap)':>12} {'|λmin|/λ+1':>11}")
print("-"*104)
for fourW in [3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]:
    M, G, N, m = weil_form(14, 44, fourW, 0.8)
    kn, kp, kz, lmin, lp1, dim = inertia(M, G)
    ratio = abs(lmin)/lp1 if lp1 and not np.isnan(lp1) else np.nan
    flag = "" if kn == 1 else "  <-- κ != 1 !!"
    print(f"{fourW:>5.1f} {m:>8} {N:>4} {dim:>4} | {kn:>7} {kp:>5} {kz:>6} | "
          f"{lmin:>+11.3e} {lp1:>+12.3e} {ratio:>11.3f}{flag}")
print("-"*104)
print("READ: κ pinned at 1 across the whole sweep (growing N, more primes) = strong direct evidence κ≤1 is")
print("ROBUST, not a low-bandwidth artifact — the single negative square persists. Watch λ_+1 (the positive")
print("floor): if it stays bounded away from 0 while λ_min stays O(1)-negative, the one negative direction is")
print("cleanly ISOLATED. If κ climbs with W, the framing needs the pole handled at each new scale. Honest either way.")
print("="*104)
