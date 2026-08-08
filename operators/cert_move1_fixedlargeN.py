#!/usr/bin/env python3
"""§839 (workflow Move 1, DECISIVE): resolve the §838 resonance-vs-stability question. PATH A found
cos(v_p,v_a)→0.9999 resonances at LOW N (5,6,7); the hypothesis is that FIXED LARGE N (out of the
resonance zone) is stable. Test: fix window LENGTH=48, 4W=3.0 (N_cap≈23), fixed N∈{12,16,20}, sweep
window HEIGHT far (bottoms 14→400). Measure cos(v_p,v_a), floor, pnorm, margin, cond(G). (cos does NOT
need real zeros — v_p,v_a are eigenvectors of the zero-independent prime/arch operators, §838 PATH D.)
GO/NO-GO: if cos stays small (~0.005) and stable across heights at fixed large N → the orthogonality
thread has legs (resonances were a low-N sampling artifact). If cos resonates/degrades even at fixed
large N → the orthogonality reduction fails, and we know it honestly."""
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
NMAX = 600; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)
RGRID = np.arange(-520.0, 520.0, 0.0025); dr = 0.0025
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

L = 48.0; fourW = 3.0; W = fourW/4; Ncap = L*2*W/np.pi
BOTTOMS = [14.0, 50.0, 100.0, 200.0, 400.0]
print("="*104)
print(f"§839 Move 1 DECISIVE — fixed LARGE N, sweep height far. length={L}, 4W={fourW}, N_cap≈{Ncap:.1f}, primes to {int(np.exp(fourW))}")
print("="*104)
for N in [12, 16, 20]:
    print("-"*104)
    print(f"N={N} (N/N_cap={N/Ncap:.2f}) FIXED; sweep window bottom 14→400:")
    print(f"  {'window':>16} {'height':>7} {'floor':>8} {'pnorm':>8} {'margin':>13} {'cos(vp,va)':>11} {'angle°':>7} {'condG':>8}")
    coss = []
    for a in BOTTOMS:
        b = a + L; CEN = list(np.linspace(a, b, N)); AP, PR, G = mats(CEN, W)
        gw = np.linalg.eigvalsh(G); condG = float(gw[-1]/max(gw[0], 1e-300))
        wf, Vf = eigh(AP, G); floor = float(wf[0]); v_a = Vf[:, 0]
        wp, Vp = eigh(PR, G); pnorm = float(wp[-1]); v_p = Vp[:, -1]
        wm, _ = eigh(AP-PR, G); margin = float(wm[0])
        cosang = abs(float(v_a @ (G @ v_p))); ang = np.degrees(np.arccos(min(1.0, cosang)))
        coss.append(cosang)
        print(f"  {f'[{a:.0f},{b:.0f}]':>16} {(a+b)/2:>7.0f} {floor:>8.4f} {pnorm:>8.4f} {margin:>+13.4e} {cosang:>11.5f} {ang:>7.2f} {condG:>8.1f}")
    print(f"  => cos range across heights: [{min(coss):.5f}, {max(coss):.5f}]  ratio {max(coss)/max(min(coss),1e-9):.1f}x  {'STABLE' if max(coss)<0.05 else 'RESONATES/DEGRADES'}")
print("="*104)
print("GO/NO-GO: all N stable (cos<0.05 across heights) → orthogonality thread has legs (resonances were low-N).")
print("          any fixed-large-N shows cos→~1 → orthogonality reduction fails; redirect honestly.")
print("="*104)
