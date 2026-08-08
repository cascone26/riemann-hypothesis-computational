#!/usr/bin/env python3
"""§840 (verification before the next workflow): confirm the two load-bearing premises of the
arch-floor-vs-prime race (§839's corrected object) before designing the workflow around them.
  P1: floor(t) grows ~ c·ln(t/2π)  (measure c across heights)
  P2: as 4W grows (more primes) at capacity-matched N, does margin stay ≥0 across heights? worst corner?
Grid: window length 30, height bottoms {14,100,400} × 4W {2.0,3.0,4.0}, N=round(0.75·N_cap), N_cap=15·4W/π.
Reports floor, pnorm, margin, cond(G). Float (verification scan)."""
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
NMAX = 800; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)
RGRID = np.arange(-520.0, 520.0, 0.004); dr = 0.004
OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)

def fej(x, W):
    y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
    out[nz] = (np.sin(y[nz])/y[nz])**2; return out
def Frow(c, r, W): return fej(r-c, W) + fej(r+c, W)

def cell(cs, W):
    n = len(cs); fourW = 4*W
    P = np.array([Frow(c, RGRID, W) for c in cs]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]
    prime = np.zeros((n, n))
    if ln.size:
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
    gw = np.linalg.eigvalsh(G); condG = float(gw[-1]/max(gw[0], 1e-300))
    floor = float(eigh(AP, G, eigvals_only=True)[0])
    pnorm = float(eigh(prime, G, eigvals_only=True)[-1])
    margin = float(eigh(AP-prime, G, eigvals_only=True)[0])
    return floor, pnorm, margin, condG

L = 30.0
print("="*100)
print("§840 VERIFY the arch-floor-vs-prime race premises (length 30, N=0.75*N_cap)")
print("="*100)
print(f"{'4W':>5} {'height':>7} {'N':>4} {'nprime':>7} {'floor':>8} {'Omega~ln':>9} {'flr/Om':>7} {'pnorm':>8} {'margin':>13} {'condG':>7}")
worst = (1e9, None)
for fourW in [2.0, 3.0, 4.0]:
    W = fourW/4; Ncap = 15*fourW/np.pi; N = max(4, round(0.75*Ncap)); npr = int((LOGN < fourW).sum())
    for a in [14.0, 100.0, 400.0]:
        b = a+L; mid = (a+b)/2; CEN = list(np.linspace(a, b, N))
        fl, pn, mg, cg = cell(CEN, W)
        Om = np.log(mid/(2*np.pi))
        if mg < worst[0]: worst = (mg, (fourW, mid, N))
        print(f"{fourW:>5.1f} {mid:>7.0f} {N:>4} {npr:>7} {fl:>8.4f} {Om:>9.3f} {fl/Om:>7.3f} {pn:>8.4f} {mg:>+13.4e} {cg:>7.1f}")
print("="*100)
print(f"P1 floor/Omega ratios above — if ~const → floor(t)≈c·ln(t/2π) confirmed.")
print(f"P2 worst-case margin = {worst[0]:+.4e} at 4W={worst[1][0]}, height={worst[1][1]:.0f}  — {'ALL ≥0 (RH-consistent)' if worst[0]>=0 else 'NEGATIVE (investigate!)'}")
print(f"   is the worst corner low-height/high-4W? {worst[1]}")
print("="*100)
