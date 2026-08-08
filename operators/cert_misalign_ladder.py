#!/usr/bin/env python3
"""§836: does the §835 orthogonality (v_prime_top ⟂ v_arch_bottom) PERSIST as window height grows?
RH ⟺ the misalignment holds uniformly. Ladder of windows at increasing height (fixed length 30,
4W=2.8 so capacity/primes fixed); at each, scan N to the razor (floor−pnorm first < 0) and report
cos(v_p,v_a) AT the razor and the margin there. If cos stays small (⟂) and margin stays 0⁺ at all
heights → orthogonality persists (RH-consistent). If cos climbs with height → the misalignment
degrades (maps the rate toward the hard side). Float."""
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

fourW = 2.80; W = fourW/4
print("="*104)
print(f"§836 ORTHOGONALITY-PERSISTENCE ladder (length 30, 4W={fourW}); does cos(v_p,v_a) stay small as height↑?")
print("="*104)
print(f"  {'window':>14} {'height':>7} {'#zeros':>7} {'razorN':>7} {'floor':>8} {'pnorm':>8} {'margin':>12} {'cos(vp,va)':>11} {'angle°':>7}")
for (a, b) in [(14.0, 44.0), (30.0, 60.0), (50.0, 80.0), (68.0, 98.0)]:
    W_ = W; nz = int(((ZEROS >= a) & (ZEROS <= b)).sum()); mid = (a+b)/2
    razorN = None; row = None
    for N in range(5, 14):
        CEN = list(np.linspace(a, b, N)); AP, PR, G = mats(CEN, W_)
        wf, Vf = eigh(AP, G); floor = float(wf[0]); v_a = Vf[:, 0]
        wp, Vp = eigh(PR, G); pnorm = float(wp[-1]); v_p = Vp[:, -1]
        wm, _ = eigh(AP-PR, G); margin = float(wm[0])
        cosang = abs(float(v_a @ (G @ v_p))); ang = np.degrees(np.arccos(min(1.0, cosang)))
        if floor - pnorm < 0 and razorN is None:
            razorN = N; row = (floor, pnorm, margin, cosang, ang); break
    if row:
        fl, pn, mg, cs_, an = row
        print(f"  {f'[{a:.0f},{b:.0f}]':>14} {mid:>7.0f} {nz:>7} {razorN:>7} {fl:>8.4f} {pn:>8.4f} {mg:>+12.4e} {cs_:>11.4f} {an:>7.2f}")
    else:
        print(f"  {f'[{a:.0f},{b:.0f}]':>14} {mid:>7.0f} {nz:>7} {'none≤13':>7}  (no razor reached within capacity)")
print("="*104)
print("READ: cos(v_p,v_a) small (angle~90°) and margin 0⁺ at all heights → orthogonality PERSISTS (RH-consistent).")
print("      cos climbing with height → misalignment degrading (a rate toward the hard side).")
print("="*104)
