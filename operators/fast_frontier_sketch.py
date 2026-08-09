#!/usr/bin/env python3
"""FAST SKETCH of the positivity frontier using float (non-rigorous) to guide the rigorous search.
This runs in ~seconds and shows us where the boundary approximately is."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
import sys

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam

NMAX = 3000; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

def analyze_margin(fourW, a=14, b=44, f_capacity=0.80, dr=0.003):
    """Fast (float) analysis of margin at 4W."""
    Ncap = 15 * fourW / np.pi
    N = max(4, round(f_capacity * Ncap))
    L = b - a
    RGRID = np.arange(-(b+60), (b+60), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)

    # Symbols
    mask = LOGN < fourW
    ln = LOGN[mask]; lw = LAMW[mask]
    Psym = np.zeros_like(RGRID)
    for u, w in zip(ln, lw): Psym += w*np.cos(u*RGRID)
    Psym *= (1/np.pi)
    S = OMEGA/(2*np.pi) - Psym

    # Weil operator
    CEN = list(np.linspace(a, b, N))
    W = fourW/4
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out

    P = np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    cosu = np.cos(np.outer(ln, RGRID))
    prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
    prime = 0.5*(prime+prime.T)

    def fc(c):
        y = W*(0.5j-c); y2 = W*(0.5j+c)
        return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in CEN])
    pole = 2.0*np.real(np.outer(fp, fp))

    Q = pole + arch - prime
    try:
        margin = float(eigh(Q, G, eigvals_only=True)[0])
    except:
        margin = -1.0

    npr = int((LOGN < fourW).sum())
    return margin, N, npr

print("FAST FRONTIER SKETCH (float, non-rigorous — for orientation only)")
print("=" * 85)
print(f"{'4W':>6} {'Ncap':>8} {'N':>4} {'N/Cap':>7} {'#primes':>7} {'float margin':>15} {'status':>20}")
print("-" * 85)

frontiers_pos = []
frontiers_neg = []

for fourW in [1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20]:
    try:
        margin, N, npr = analyze_margin(fourW)
        Ncap = 15*fourW/np.pi
        status = "✓ POSITIVE" if margin > 1e-6 else ("THIN" if margin > -1e-4 else "✗ NEGATIVE")
        print(f"{fourW:>6.2f} {Ncap:>8.2f} {N:>4} {N/Ncap:>7.2f} {npr:>7} {margin:>+15.8e}  {status:>20}")
        if margin > 1e-6:
            frontiers_pos.append(fourW)
        else:
            frontiers_neg.append(fourW)
    except Exception as e:
        print(f"{fourW:>6.2f} {'?':>8} {'?':>4} {'?':>7} {'?':>7} {'ERROR':>15}  {str(e)[:20]:>20}")

print("=" * 85)
if frontiers_pos:
    print(f"\nLast POSITIVE (float): 4W ≈ {max(frontiers_pos):.2f}")
if frontiers_neg:
    print(f"First NEGATIVE (float): 4W ≈ {min(frontiers_neg):.2f}")

if frontiers_pos and frontiers_neg:
    interval_lo = max(frontiers_pos)
    interval_hi = min(frontiers_neg)
    print(f"\nRigorous bisection should probe: [{interval_lo:.2f}, {interval_hi:.2f}]")
    print(f"Expected frontier (approximate): 4W ≈ {(interval_lo + interval_hi)/2:.2f}")
else:
    print(f"\nAll positive: frontier is beyond 4W ≈ {max(frontiers_pos):.2f}")

print("\n" + "=" * 85)
print("This is a FLOAT sketch only. Use rigid arb bisection for certified frontier.")
