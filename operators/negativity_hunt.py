#!/usr/bin/env python3
"""§852: broad NEGATIVITY HUNT. Scan a wide grid of (window height, 4W, N); compute the full Weil margin
= λmin(pole+arch−prime, G) and G-conditioning. Report the MOST NEGATIVE margin found and whether ANY
window has a robustly-negative margin (margin < −1e-6 AND cond(G) healthy < 100 → not a conditioning
artifact). If RH holds, none should be robustly negative. A robust negative = counterexample or bug (BIG).
Everywhere ≥0 = strong RH-consistent evidence across a wide sweep. Float."""
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

def margin_cond(a, b, W, N, dr=0.005):
    fourW = 4*W
    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    CEN = list(np.linspace(a, b, N))
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    Pm = np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G = (Pm @ Pm.T)*dr
    gw = np.linalg.eigvalsh(G); cond = float(gw[-1]/max(gw[0], 1e-300))
    arch = (Pm*OMEGA[None, :]) @ Pm.T * dr/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]
    cosu = np.cos(np.outer(ln, RGRID)); prime = np.zeros((N, N))
    for j in range(N):
        gp = ((Pm[j][None, :]*Pm) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
    prime = 0.5*(prime+prime.T)
    def fc(c):
        y = W*(0.5j-c); y2 = W*(0.5j+c); return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in CEN]); pole = 2.0*np.real(np.outer(fp, fp))
    return float(eigh(pole+arch-prime, G, eigvals_only=True)[0]), cond

print("="*90)
print("§852 NEGATIVITY HUNT — does the Weil margin EVER go robustly negative (cond-healthy)?")
print("="*90)
robust_neg = []
worst = (1e9, None); count = 0
for a in [14, 30, 60, 120, 250, 500]:
    for L in [20, 30, 40]:
        b = a+L
        for fourW in [1.6, 2.0, 2.5, 3.0, 3.5, 4.0]:
            W = fourW/4; Ncap = L*2*W/np.pi
            for frac in [0.5, 0.7, 0.9, 1.05]:
                N = max(4, round(frac*Ncap))
                if N > 24: continue
                mg, cond = margin_cond(a, b, W, N); count += 1
                if mg < worst[0]: worst = (mg, (a, b, fourW, N, cond))
                if mg < -1e-6 and cond < 100:
                    robust_neg.append((a, b, fourW, N, mg, cond))
print(f"scanned {count} (window,4W,N) configs.")
print(f"MOST NEGATIVE margin: {worst[0]:+.4e} at window=[{worst[1][0]},{worst[1][1]}], 4W={worst[1][2]}, N={worst[1][3]}, cond={worst[1][4]:.1f}")
if robust_neg:
    print(f"\n*** {len(robust_neg)} ROBUSTLY-NEGATIVE (margin<-1e-6, cond<100) configs found — INVESTIGATE ***")
    for r in robust_neg[:10]: print("   ", r)
else:
    print(f"\nNO robustly-negative margin found (all margin ≥ −1e-6 OR only-negative-when-ill-conditioned).")
    print("Strong RH-consistent evidence: Weil positivity holds across the entire scanned grid.")
print("="*90)
