#!/usr/bin/env python3
"""§966 DECISIVE: does #pos(arch−prime) track 2·RvM(window) as the window MOVES and WIDENS?  (follows §965)

§965: at window [14,44], saturated #pos = 15 ≈ 2·RvM = 14.14, and #neg≡1. If the positive rank IS the
zero-count, then across DIFFERENT windows (different center T and width) #pos must track 2·RvM(window) each
time, while #neg stays 1. That nails 'the positive spectrum counts the zeros in the window' as a genuine
Weyl identity, not a coincidence at one window. Run each window at a SATURATED bandwidth (4W=7, saturation
reached by 4W=5 per §964/965); report #neg, #pos, and 2·RvM. HANDLE (own code, honest)."""
import numpy as np
from scipy.special import digamma

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam
NMAX = 8000; LAM = vonmangoldt(NMAX)
NS = np.arange(2, NMAX+1); LOGN = np.log(NS); LAMv = LAM[2:]

def RvM(T):
    if T <= 0: return 0.0
    return (T/(2*np.pi))*np.log(T/(2*np.pi)) - T/(2*np.pi) + 7/8

def build(cen_lo, cen_hi, fourW, frac, dr=0.004):
    W = fourW/4; L = cen_hi-cen_lo; N = max(4, round(frac*L*2*W/np.pi))
    RG = np.arange(-(cen_hi+120), (cen_hi+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RG)) - np.log(np.pi)
    CEN = np.linspace(cen_lo, cen_hi, N)
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    P = np.array([fej(RG-c)+fej(RG+c) for c in CEN]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; wts = LAMv[mask]/np.sqrt(NS[mask])
    cosu = np.cos(np.outer(ln, RG)); prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*wts[None, :]).sum(axis=1)
    return arch, 0.5*(prime+prime.T), G, N

def inertia(M, G):
    gw, gv = np.linalg.eigh(G); keep = gw > gw.max()*1e-12
    B = gv[:, keep]/np.sqrt(gw[keep]); ev = np.sort(np.linalg.eigvalsh(B.T @ M @ B))
    s = max(1.0, np.abs(ev).max()); tol = 1e-9*s
    return int((ev < -tol).sum()), int((ev > tol).sum()), int((np.abs(ev) <= tol).sum())

windows = [(14,30), (14,44), (14,64), (30,60), (20,80), (50,90), (80,120), (100,150)]
print("="*104)
print("§966 does #pos(arch−prime) track 2·RvM(window)?  saturated bandwidth 4W=7.  ζ (q=1).")
print("="*104)
print(f"{'window':>12} {'width':>6} {'N':>4} | {'#neg':>5} {'#pos':>5} {'#zero':>6} | "
      f"{'RvM(win)':>9} {'2·RvM':>7} {'#pos−2RvM':>10}")
print("-"*104)
for (lo, hi) in windows:
    arch, prime, G, N = build(lo, hi, 7.0, 0.8)
    dn, dp, dz = inertia(arch-prime, G)
    z = RvM(hi) - RvM(lo)
    print(f"[{lo:>3},{hi:>3}]{'':>3} {hi-lo:>6} {N:>4} | {dn:>5} {dp:>5} {dz:>6} | "
          f"{z:>9.3f} {2*z:>7.2f} {dp-2*z:>+10.2f}")
print("-"*104)
print("READ: if #pos ≈ 2·RvM(window) across ALL windows (tracking as the window moves higher and widens) with")
print("#neg≡1, then the positive rank of arch−prime IS the Weyl zero-count of the window — 'zeros = positive")
print("eigenvalues' as a genuine identity, and κ=1 is the single-unit overshoot (the pole). If #pos drifts away")
print("from 2·RvM, the §965 match was a one-window coincidence. Honest either way — report the residual column.")
print("="*104)
