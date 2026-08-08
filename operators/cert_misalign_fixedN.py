#!/usr/bin/env python3
"""§837: de-confound §836. Measure cos(v_p,v_a) at FIXED N (⟹ fixed capacity-fraction, since all windows
are length 30 at 4W=2.8 → identical N_cap≈13.4) across a height ladder. If cos still climbs with height
at fixed N → the misalignment degradation is a REAL height effect (not razorN drifting toward capacity).
If cos flattens at fixed N → §836's climb was the capacity-proximity confound. Float."""
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
RGRID = np.arange(-320.0, 320.0, 0.0025); dr = 0.0025
OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)

def fej(x, W):
    y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
    out[nz] = (np.sin(y[nz])/y[nz])**2; return out
def Frow(c, r, W): return fej(r-c, W) + fej(r+c, W)

def mats(cs, W):
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
    return pole+arch, prime, G

fourW = 2.80; W = fourW/4; Ncap = 30*2*W/np.pi
WINDOWS = [(14.0, 44.0), (30.0, 60.0), (50.0, 80.0), (68.0, 98.0)]
print("="*104)
print(f"§837 FIXED-N misalignment de-confound (length 30, 4W={fourW}, N_cap≈{Ncap:.1f}); cos vs height at fixed N")
print("="*104)
for N in [8, 10, 12]:
    print("-"*104)
    print(f"N={N} (N/N_cap={N/Ncap:.2f}) held FIXED across heights:")
    print(f"  {'window':>14} {'height':>7} {'floor':>8} {'pnorm':>8} {'floor-pnorm':>12} {'margin':>12} {'cos(vp,va)':>11} {'angle°':>7} {'condG':>8}")
    for (a, b) in WINDOWS:
        CEN = list(np.linspace(a, b, N)); AP, PR, G = mats(CEN, W)
        gw = np.linalg.eigvalsh(G); condG = float(gw[-1]/max(gw[0], 1e-300))
        wf, Vf = eigh(AP, G); floor = float(wf[0]); v_a = Vf[:, 0]
        wp, Vp = eigh(PR, G); pnorm = float(wp[-1]); v_p = Vp[:, -1]
        wm, _ = eigh(AP-PR, G); margin = float(wm[0])
        cosang = abs(float(v_a @ (G @ v_p))); ang = np.degrees(np.arccos(min(1.0, cosang)))
        print(f"  {f'[{a:.0f},{b:.0f}]':>14} {(a+b)/2:>7.0f} {floor:>8.4f} {pnorm:>8.4f} {floor-pnorm:>12.5f} {margin:>+12.4e} {cosang:>11.4f} {ang:>7.2f} {condG:>8.1f}")
print("="*104)
print("READ: at FIXED N, if cos still climbs with height → real height effect on misalignment.")
print("      if cos ~flat across heights → §836's climb was razorN→capacity confound (the honest resolution).")
print("="*104)
