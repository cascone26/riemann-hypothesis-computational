#!/usr/bin/env python3
"""§968 THE RESOLUTION-THRESHOLD LAW 4W_onset(T) — quantify the prime–zero duality behind κ=1.  (HANDLE)

§967: κ jumps 0→1 at a bandwidth ONSET that grows with window height T. Quantify it: for narrow windows
centered at increasing T, sweep 4W on a fine grid and record the smallest 4W where κ=1 (and holds). Then
compare 4W_onset(T) to candidate laws — LINEAR (4W~cT) vs LOG (4W~c·lnT) vs the explicit-formula pairing.
The Weil prime sum uses primes with log n < 4W; the archimedean phase over height T advances like
∫Ω ~ T·ln(T/2π). If 4W_onset tracks a clean function of T, that is a real quantitative bridge: 'to see the
single negative square (pole) at height T you need primes up to e^{4W_onset(T)}'. HANDLE (own code, honest)."""
import numpy as np
from scipy.special import digamma

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam
NMAX = 4000; LAM = vonmangoldt(NMAX)
NS = np.arange(2, NMAX+1); LOGN = np.log(NS); LAMv = LAM[2:]

def kap_at(cen_lo, cen_hi, fourW, frac=0.8, dr=0.005):
    W = fourW/4; L = cen_hi-cen_lo; N = max(4, round(frac*L*2*W/np.pi))
    RG = np.arange(-(cen_hi+40), (cen_hi+40), dr)
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
    M = arch - 0.5*(prime+prime.T)
    gw, gv = np.linalg.eigh(G); keep = gw > gw.max()*1e-12
    B = gv[:, keep]/np.sqrt(gw[keep]); ev = np.linalg.eigvalsh(B.T @ M @ B)
    s = max(1.0, np.abs(ev).max()); return int((ev < -1e-9*s).sum())

grid = np.arange(3.0, 8.01, 0.5)
centers = [20, 30, 40, 50, 60]
half = 9  # narrow window half-width
print("="*92)
print("§968 onset bandwidth 4W_onset(T) where κ jumps to 1. narrow windows [T-9,T+9].")
print("="*92)
print(f"{'T':>4} | " + "".join(f"{w:>4.1f}" for w in grid) + " | 4W_onset")
print("-"*92)
onsets = []
for T in centers:
    row = []; onset = None
    for w in grid:
        k = kap_at(T-half, T+half, w)
        row.append(k)
        if onset is None and k >= 1: onset = w
    onsets.append((T, onset))
    os = f"{onset:.1f}" if onset else ">8"
    print(f"{T:>4} | " + "".join(f"{k:>4}" for k in row) + f" | {os}")
print("-"*92)
print("candidate laws (fit onset vs T):")
Ts = np.array([t for t, o in onsets if o]); Os = np.array([o for t, o in onsets if o], float)
if len(Ts) >= 2:
    # linear 4W = a*T + b
    A = np.polyfit(Ts, Os, 1); pred_lin = np.polyval(A, Ts)
    # log 4W = a*lnT + b
    B = np.polyfit(np.log(Ts), Os, 1); pred_log = np.polyval(B, np.log(Ts))
    r_lin = np.sqrt(np.mean((Os-pred_lin)**2)); r_log = np.sqrt(np.mean((Os-pred_log)**2))
    print(f"  LINEAR  4W_onset ≈ {A[0]:.4f}·T + {A[1]:.3f}   rms={r_lin:.3f}")
    print(f"  LOG     4W_onset ≈ {B[0]:.4f}·lnT + {B[1]:.3f}  rms={r_log:.3f}")
    print(f"  ratio T/4W_onset per point: " + ", ".join(f"{t}:{t/o:.1f}" for t, o in zip(Ts, Os)))
    print(f"  => {'LINEAR' if r_lin<r_log else 'LOG'} fits better. explicit-formula reading:")
    print(f"     primes needed up to e^(4W_onset); zeros to height T need log-prime reach ~ {'T (linear)' if r_lin<r_log else 'ln T (log)'}.")
print("READ: a clean law connects κ=1's resolution threshold to the prime–zero explicit-formula duality —")
print("quantifying WHY primes to a certain height are exactly what's needed to expose the single pole-negative.")
print("="*92)
