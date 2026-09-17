#!/usr/bin/env python3
"""§967 κ≤1 UNIVERSALITY AT SATURATION — the CORRECT test after §966's kill.  (HANDLE)

§966 showed κ at a FIXED bandwidth is confounded (high-T windows under-saturate, #neg spuriously 0). The real
statement is: κ=1 is a SATURATION phenomenon — once 4W is large enough for the window's height, does EVERY
window (all along the critical line) land at κ=1? Sweep 4W upward per window; print the κ trajectory; a window
whose κ has stopped changing (and equals its neighbors' saturated value) is saturated. RH-content: κ_sat=1 for
every window. A window saturating at κ=2 would be a real anomaly (or bug) worth chasing. HANDLE (own code)."""
import numpy as np
from scipy.special import digamma

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam
NMAX = 8000; LAM = vonmangoldt(NMAX)
NS = np.arange(2, NMAX+1); LOGN = np.log(NS); LAMv = LAM[2:]

def build(cen_lo, cen_hi, fourW, frac, dr=0.004):
    W = fourW/4; L = cen_hi-cen_lo; N = max(4, round(frac*L*2*W/np.pi))
    RG = np.arange(-(cen_hi+60), (cen_hi+60), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RG)) - np.log(np.pi)
    CEN = np.linspace(cen_lo, cen_hi, N)
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    P = np.array([fej(RG-c)+fej(RG+c) for c in CEN]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; wts = LAMv[mask]/np.sqrt(NS[mask])
    cosu = np.cos(np.outer(ln, RG)); prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*wts[None, :]).sum(axis=1)
    return arch, 0.5*(prime+prime.T), G, N

def kap(M, G):
    gw, gv = np.linalg.eigh(G); keep = gw > gw.max()*1e-12
    B = gv[:, keep]/np.sqrt(gw[keep]); ev = np.linalg.eigvalsh(B.T @ M @ B)
    s = max(1.0, np.abs(ev).max()); return int((ev < -1e-9*s).sum())

windows = [(14,44), (30,60), (45,75), (60,90)]
sweeps   = [4.0, 5.0, 6.0, 7.0, 8.0]
print("="*96)
print("§967 κ(arch−prime) trajectory vs bandwidth, per window — does each SATURATE at κ=1 along the line?")
print("="*96)
header = "     window |" + "".join(f"  4W={w:.0f}" for w in sweeps) + "  | saturated κ"
print(header); print("-"*96)
for (lo, hi) in windows:
    traj = []
    for fourW in sweeps:
        arch, prime, G, N = build(lo, hi, fourW, 0.8)
        traj.append(kap(arch-prime, G))
    # saturated = last value if last two agree, else 'not sat (need higher 4W)'
    sat = str(traj[-1]) if traj[-1] == traj[-2] else f"~{traj[-1]}?(rising)"
    row = f"  [{lo:>3},{hi:>3}] |" + "".join(f"  {t:>4}" for t in traj) + f"  | {sat}"
    print(row)
print("-"*96)
print("READ: each window's κ should RISE from 0 (under-resolved) to 1 as 4W grows, then HOLD at 1. If every")
print("window along the line saturates at κ=1 (last two bandwidths agree at 1), that is κ≤1 universality across")
print("the critical strip at saturation — the robust content behind RH in this framing. A window locking at κ=2")
print("would be a genuine anomaly. Honest: report the trajectories; 'still rising' = needs higher 4W, not κ=that.")
print("="*96)
