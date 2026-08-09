#!/usr/bin/env python3
"""§851: solidify the §850 law (#near-null Weil eigenvalues = max(0, N − #zeros-in-window)) across MANY
windows with REAL zero counts (mpmath.zetazero). Vary window height and N; count zeros in each window
and the near-null Weil eigenvalues; check the law holds broadly."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
import mpmath as mp

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam
NMAX = 3000; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

# precompute zero heights up to ~600 (enough for windows tested)
mp.mp.dps = 20
ZTAB = []
n = 1
while True:
    g = float(mp.im(mp.zetazero(n)))
    ZTAB.append(g)
    if g > 640: break
    n += 1
ZTAB = np.array(ZTAB)

def zeros_in(a, b): return int(((ZTAB >= a) & (ZTAB <= b)).sum())

def near_null(a, b, W, thresh, dr=0.004):
    fourW = 4*W; L = b-a
    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    Nn = max(4, round(0.9*L*2*W/np.pi))
    CEN = list(np.linspace(a, b, Nn))
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    Pm = np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G = (Pm @ Pm.T)*dr
    arch = (Pm*OMEGA[None, :]) @ Pm.T * dr/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]
    cosu = np.cos(np.outer(ln, RGRID)); prime = np.zeros((Nn, Nn))
    for j in range(Nn):
        gp = ((Pm[j][None, :]*Pm) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
    prime = 0.5*(prime+prime.T)
    def fc(c):
        y = W*(0.5j-c); y2 = W*(0.5j+c); return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in CEN]); pole = 2.0*np.real(np.outer(fp, fp))
    w = eigh(pole+arch-prime, G, eigvals_only=True)
    return Nn, int((w < thresh).sum()), float(w[0])

print("="*96)
print("§851 verify #near-null = max(0, N − #zeros-in-window) across windows (real zeros via mpmath)")
print("="*96)
print(f"{'window':>14} {'#zeros':>7} {'4W':>4} {'N':>3} {'#near-null':>11} {'N−#zeros':>9} {'match?':>7} {'margin':>11}")
ok_all = True
for (a, b, fourW) in [(14,44,3.0),(14,44,4.0),(50,80,3.0),(50,80,4.0),(150,180,3.0),(150,180,4.0),(300,330,4.0)]:
    W = fourW/4; nz = zeros_in(a, b)
    thr = max(1e-3, 0.02)   # near-null threshold (well below the O(0.1-1) bulk)
    Nn, nn, mg = near_null(a, b, W, thr)
    pred = max(0, Nn - nz); match = (nn == pred)
    ok_all = ok_all and match
    print(f"{f'[{a},{b}]':>14} {nz:>7} {fourW:>4.1f} {Nn:>3} {nn:>11} {pred:>9} {str(match):>7} {mg:>+11.3e}", flush=True)
print("="*96)
print(f"LAW #near-null = max(0,N−#zeros): {'HOLDS across all tested windows' if ok_all else 'FAILS somewhere (see above)'}")
print("="*96)
