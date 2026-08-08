#!/usr/bin/env python3
"""§833: THE decisive test of the §832 interpretation — is the finite Weil razor set by the WINDOW'S
ZERO-COUNT? Hold window LENGTH fixed (30) and 4W fixed (2.8 → identical primes, identical frame
capacity N_cap≈13.4), vary only the window HEIGHT so ONLY the zero-count changes:
   [14,44]=8 zeros, [40,70]=11 zeros, [68,98]=12 zeros.
Sweep N; the margin-collapse onset (first N with margin<1e-2) should land at N≈#zeros if the razor is
zero-count-set. Length/primes/capacity are identical across the three, so a moving onset PROVES it's the
zero-count, not length or conditioning. Float (fast); the §832 arb ladder already certified the values."""
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
NMAX = 400; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)
RGRID = np.arange(-260.0, 260.0, 0.0025); dr = 0.0025
OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
ZEROS = np.array([
    14.134725142, 21.022039639, 25.010857580, 30.424876126, 32.935061588,
    37.586178159, 40.918719012, 43.327073281, 48.005150881, 49.773832478,
    52.970321478, 56.446247697, 59.347044003, 60.831778525, 65.112544048,
    67.079810529, 69.546401711, 72.067157674, 75.704690699, 77.144840069,
    79.337375020, 82.910380854, 84.735492981, 87.425274613, 88.809111208,
    92.491899271, 94.651344041, 95.870634228, 98.831194218])

def fej(x, W):
    y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
    out[nz] = (np.sin(y[nz])/y[nz])**2; return out
def Frow(c, r, W): return fej(r-c, W) + fej(r+c, W)

def build(cs, W):
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
    return pole+arch-prime, G

def diag(cs, W):
    Q, G = build(cs, W)
    gw = np.linalg.eigvalsh(G); lg = float(gw[0]); cg = float(gw[-1]/max(gw[0], 1e-300))
    return float(eigh(Q, G, eigvals_only=True)[0]), lg, cg

fourW = 2.80; W = fourW/4; Ncap = 30*2*W/np.pi
print("="*100)
print(f"§833 ZERO-COUNT RAZOR TEST — length=30, 4W={fourW} (N_cap≈{Ncap:.1f}, primes to {int(np.exp(fourW))}) FIXED;")
print("only the window HEIGHT (⟹ zero-count) varies. Onset (margin<1e-2) should track N=#zeros.")
print("="*100)
windows = [(14.0, 44.0), (40.0, 70.0), (68.0, 98.0)]
for (a, b) in windows:
    nz = int(((ZEROS >= a) & (ZEROS <= b)).sum())
    print("-"*100)
    print(f"window [{a:.0f},{b:.0f}]  #zeros={nz}   (predicted collapse onset: N≈{nz})")
    print(f"  {'N':>4} {'margin':>14} {'λmin(G)':>10} {'cond(G)':>9}  {'health':>10}")
    onset = None
    for N in range(4, min(int(Ncap)+1, nz+4)+1):
        CEN = list(np.linspace(a, b, N)); m, lg, cg = diag(CEN, W)
        health = "G-healthy" if (lg > 0.3 and cg < 40) else ("G-loading" if cg < 1e3 else "G-DEGEN")
        flag = ""
        if onset is None and m < 1e-2:
            onset = N; flag = "  <-- ONSET"
        print(f"  {N:>4} {m:>+14.5e} {lg:>10.3f} {cg:>9.1f}  {health:>10}{flag}")
    print(f"  => collapse onset at N={onset}  vs  #zeros={nz}  :  {'MATCH' if onset==nz else ('off by '+str((onset or -1)-nz))}")
print("="*100)
print("If onset tracks #zeros across the 3 windows (8/11/12) at FIXED length+4W+primes → the finite Weil")
print("razor is PROVABLY set by the window's zero-count — the tightest finite RH-shadow (zeros-on-line ⟺")
print("Weil form marginally positive at critical sampling).")
print("="*100)
