#!/usr/bin/env python3
"""§855: do the RH-critical NEAR-NULL eigenvalues have arithmetic content — do they track the window's
ZERO SPACINGS? The near-null modes (§850) are band-limited functions in the span that best vanish at the
window's zeros; their Weil value = 2Σ_γ φ̂(γ)² for the residual, which interpolation theory ties to the
zeros' configuration. Test: at N = #zeros + 1 (exactly ONE near-null mode = the margin, kept below
G-capacity so no conditioning confound), correlate that eigenvalue with the window's min & mean zero
spacing. Real zeros via mpmath. Full Weil form (no §853 confound). Honest if no clean law."""
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
NMAX = 4000; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

mp.mp.dps = 20
ZT = []
k = 1
while True:
    g = float(mp.im(mp.zetazero(k))); ZT.append(g)
    if g > 360: break
    k += 1
ZT = np.array(ZT)

def zeros_in(a, b): return ZT[(ZT >= a) & (ZT <= b)]

def weil_spectrum(a, b, W, N, dr=0.004):
    fourW = 4*W
    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    CEN = list(np.linspace(a, b, N))
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    P = np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G = (P @ P.T)*dr
    gw = np.linalg.eigvalsh(G); cond = float(gw[-1]/max(gw[0], 1e-300))
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]
    cosu = np.cos(np.outer(ln, RGRID)); prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
    prime = 0.5*(prime+prime.T)
    def fc(c):
        y = W*(0.5j-c); y2 = W*(0.5j+c); return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in CEN]); pole = 2.0*np.real(np.outer(fp, fp))
    w = eigh(pole+arch-prime, G, eigvals_only=True)
    return w, cond

print("="*104)
print("§855 near-null eigenvalue vs window zero-spacings — arithmetic content? (N=#zeros+1, 1 near-null mode)")
print("="*104)
print(f"{'window':>12} {'4W':>4} {'#z':>3} {'N':>3} {'cond':>6} {'margin(near-null)':>17} {'min_gap':>8} {'mean_gap':>8} {'m/min_gap^?':>12}")
rows = []
for (a, b, fourW) in [(14,44,4.0),(30,60,4.0),(60,90,4.0),(100,130,5.0),(150,180,5.0)]:
    W = fourW/4; zs = zeros_in(a, b); nz = len(zs)
    if nz < 2: continue
    gaps = np.diff(np.sort(zs)); min_gap = float(gaps.min()); mean_gap = float(gaps.mean())
    N = nz + 1; Ncap = (b-a)*2*W/np.pi
    if N > Ncap - 0.5:  # keep below capacity
        continue
    w, cond = weil_spectrum(a, b, W, N)
    margin = float(w[0])
    rows.append((a, b, nz, N, cond, margin, min_gap, mean_gap))
    print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {nz:>3} {N:>3} {cond:>6.1f} {margin:>+17.4e} {min_gap:>8.3f} {mean_gap:>8.3f} {margin/min_gap:>12.3e}")
print("="*104)
# look for a clean scaling: margin vs min_gap (log-log)
if len(rows) >= 3:
    mg = np.array([r[5] for r in rows]); mn = np.array([r[6] for r in rows]); mnn = np.array([r[7] for r in rows])
    pos = mg > 0
    if pos.sum() >= 3:
        for name, x in [("min_gap", mn[pos]), ("mean_gap", mnn[pos])]:
            b1, b0 = np.polyfit(np.log(x), np.log(mg[pos]), 1)
            r2 = 1 - np.sum((np.log(mg[pos]) - (b1*np.log(x)+b0))**2)/np.sum((np.log(mg[pos])-np.log(mg[pos]).mean())**2)
            print(f"log(margin) ~ {b1:.2f}·log({name}) + {b0:.2f}   R²={r2:.3f}")
print("If margin ~ (min_gap)^p cleanly (high R²) → the RH-critical eigenvalue has an arithmetic law in the zero spacings.")
print("="*104)
