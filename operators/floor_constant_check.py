#!/usr/bin/env python3
"""§843: is the floor coefficient the EXACT constant 1/(2π)=0.159155? PATH C measured floor/Ω→0.159.
The arch term carries a 1/2π, so the floor (min arch eigenvalue) should asymptote to (1/2π)·(effective
Ω over the window). Test: compute floor(t) at high t, ratio to (1/2π)·ln(t/2π). If →1, the floor law's
leading coefficient is EXACTLY 1/(2π) — an exact grounding, not a fit."""
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
NMAX = 5000; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

def floor_and_Omegabar(a, b, W, dr=0.004):
    L = b-a; RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    N = max(4, round(0.75*L*2*W/np.pi)); CEN = list(np.linspace(a, b, N))
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    P = np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    def fc(c):
        y = W*(0.5j-c); y2 = W*(0.5j+c); return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in CEN]); pole = 2.0*np.real(np.outer(fp, fp))
    floor = float(eigh(pole+arch, G, eigvals_only=True)[0])
    # effective Omega over window (mean of Omega on [a,b])
    win = (RGRID >= a) & (RGRID <= b)
    Omega_bar = float(OMEGA[win].mean())
    return floor, Omega_bar, N

inv2pi = 1/(2*np.pi)
print("="*94)
print(f"§843 is the floor coefficient EXACTLY 1/(2π)={inv2pi:.6f}?  (length 40, 4W=3)")
print("="*94)
print(f"{'height t':>9} {'N':>4} {'floor':>10} {'ln(t/2pi)':>10} {'floor/ln':>10} {'(1/2pi)*Ombar':>14} {'floor/[(1/2pi)ln]':>18}")
for a in [60.0, 200.0, 500.0, 1500.0, 4000.0]:
    b = a+40; mid = (a+b)/2; W = 0.75
    fl, Omb, N = floor_and_Omegabar(a, b, W)
    lnt = np.log(mid/(2*np.pi))
    print(f"{mid:>9.0f} {N:>4} {fl:>10.5f} {lnt:>10.4f} {fl/lnt:>10.5f} {inv2pi*Omb:>14.5f} {fl/(inv2pi*lnt):>18.5f}", flush=True)
print("="*94)
print(f"if floor/[(1/2pi)ln(t/2pi)] → 1  →  floor(t) ~ (1/2π)·ln(t/2π) EXACTLY (coefficient grounded, not fitted).")
print(f"and floor ≈ (1/2π)·Ω̄(window)  →  the floor IS the window-averaged arch density /2π.")
print("="*94)
