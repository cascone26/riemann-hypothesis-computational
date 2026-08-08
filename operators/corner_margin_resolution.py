#!/usr/bin/env python3
"""§841: is the corner margin (1.2e-9 at [14,44], 4W=4) a REAL tiny positive margin asymptoting to 0+
(RH-consistent) or float noise around 0? Recompute at 3 integration resolutions (dr=0.004/0.002/0.001).
Stable across resolutions → real; jumps → noise/undetermined. Complements the workflow's PATH A."""
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
NMAX = 2000; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

def make(dr):
    RGRID = np.arange(-300.0, 300.0, dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    return RGRID, OMEGA

def fej(x, W):
    y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
    out[nz] = (np.sin(y[nz])/y[nz])**2; return out
def Frow(c, r, W): return fej(r-c, W) + fej(r+c, W)

def margin(cs, W, dr):
    RGRID, OMEGA = make(dr); n = len(cs); fourW = 4*W
    P = np.array([Frow(c, RGRID, W) for c in cs]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]
    prime = np.zeros((n, n))
    cosu = np.cos(np.outer(ln, RGRID))
    for j in range(n):
        gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
    prime = 0.5*(prime+prime.T)
    def fc(c):
        y = W*(0.5j-c); y2 = W*(0.5j+c)
        return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in cs]); pole = 2.0*np.real(np.outer(fp, fp))
    AP = pole+arch
    fl = float(eigh(AP, G, eigvals_only=True)[0])
    pn = float(eigh(prime, G, eigvals_only=True)[-1])
    mg = float(eigh(AP-prime, G, eigvals_only=True)[0])
    return fl, pn, mg

print("="*80)
print("§841 corner margin resolution-stability: [14,44], 4W=4.0, N=14")
print("real tiny +margin (stable across dr) or float noise (jumps)?")
print("="*80)
W = 1.0; CEN = list(np.linspace(14, 44, 14))
print(f"{'dr':>8} {'floor':>10} {'pnorm':>10} {'margin':>16}")
for dr in [0.004, 0.002, 0.001]:
    fl, pn, mg = margin(CEN, W, dr)
    print(f"{dr:>8.4f} {fl:>10.5f} {pn:>10.5f} {mg:>+16.6e}", flush=True)
print("="*80)
print("stable across dr (same sign+order) → real margin; wild swings → float noise / undetermined")
print("="*80)
