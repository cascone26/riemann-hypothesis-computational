#!/usr/bin/env python3
"""PATH C — NEGATIVE-SET CONCENTRATION: Refine §846. The lift is band-limited functions
AVOIDING the symbol's negative set {r∈window : S(r)<0}.

This script computes:
  1. measure/structure of {S<0} within the window (how much is negative, how many components)
  2. top eigenvalue c_neg of the concentration operator onto {S<0}
     (band-limited functions restricted to the negative set, in G-metric)
  3. does lift correlate with (1−c_neg) or with negative-set measure?

HYPOTHESIS: If band-limited functions of type 2·4W CAN'T concentrate on {S<0} (c_neg
bounded away from 1), that's WHY the lift compensates — an uncertainty-principle mechanism.
"""
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

def analyze_negative_set(a, b, W, dr=0.004):
    """Compute negative-set structure and concentration."""
    fourW = 4*W; L = b-a
    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)

    # Build Weil symbol S(r) = (1/2π)Ω(r) - P(r)
    arch_sym = OMEGA/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]
    P = np.zeros_like(RGRID)
    for u, w in zip(ln, lw):
        P += w*np.cos(u*RGRID)
    P *= (1.0/np.pi)
    S = arch_sym - P

    # Window and negative set
    win_idx = (RGRID >= a) & (RGRID <= b)
    S_win = S[win_idx]
    RGRID_win = RGRID[win_idx]

    # Negative set structure
    neg_idx = S_win < 0
    S_neg = S_win[neg_idx]
    RGRID_neg = RGRID_win[neg_idx]

    # Negative-set measures
    neg_count = neg_idx.sum()
    neg_measure = float(neg_count * dr)  # total length of negative set
    neg_measure_frac = float(neg_count / len(S_win))  # fraction of window
    S_min = float(S_win.min())
    S_neg_L1 = float(np.maximum(0, -S_win).sum() * dr)  # L1 mass of negative part

    # Find connected components of negative set (within window)
    components = []
    if neg_count > 0:
        diffs = np.diff(np.where(neg_idx)[0])
        breaks = np.where(diffs > 1)[0] + 1
        splits = np.split(np.where(neg_idx)[0], breaks)
        for split in splits:
            if len(split) > 0:
                r_min = RGRID_win[split[0]]
                r_max = RGRID_win[split[-1]]
                components.append((r_min, r_max, len(split)*dr))
    num_components = len(components)

    # Build basis functions and Gram matrix
    N = max(4, round(0.75*L*2*W/np.pi)); CEN = list(np.linspace(a, b, N))
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out

    Pm = np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN])
    G = (Pm @ Pm.T)*dr

    # Full Weil margin (for reference)
    arch = (Pm*OMEGA[None, :]) @ Pm.T * dr/(2*np.pi)
    cosu = np.cos(np.outer(ln, RGRID)); prime = np.zeros((N, N))
    for j in range(N):
        gp = ((Pm[j][None, :]*Pm) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
    prime = 0.5*(prime+prime.T)
    def fc(c):
        y = W*(0.5j-c); y2 = W*(0.5j+c); return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in CEN]); pole = 2.0*np.real(np.outer(fp, fp))
    margin = float(eigh(pole+arch-prime, G, eigvals_only=True)[0])

    # Full-window concentration (prolate: top eig of A_window = 2 int_window F_j F_k)
    A_window = (Pm[:, win_idx] @ Pm[:, win_idx].T) * dr * 2.0
    conc_window = float(eigh(A_window, G, eigvals_only=True)[-1])

    # Negative-set concentration: A_neg = 2 int_{S<0} F_j F_k
    if neg_count > 0:
        A_neg = (Pm[:, win_idx][:, neg_idx] @ Pm[:, win_idx][:, neg_idx].T) * dr * 2.0
        c_neg = float(eigh(A_neg, G, eigvals_only=True)[-1])
    else:
        c_neg = 0.0

    # Lift (margin - S_min)
    lift = margin - S_min

    return {
        'S_min': S_min,
        'margin': margin,
        'lift': lift,
        'neg_measure': neg_measure,
        'neg_measure_frac': neg_measure_frac,
        'num_components': num_components,
        'S_neg_L1': S_neg_L1,
        'conc_window': conc_window,
        'c_neg': c_neg,
        'N': N,
        'components': components
    }

print("="*130)
print("PATH C — NEGATIVE-SET CONCENTRATION ANALYSIS")
print("="*130)
print(f"{'window':>12} {'4W':>4} {'S_min':>9} {'margin':>12} {'lift':>9} {'neg_frac':>9} {'#comp':>5} {'Sneg_L1':>9}")
print(f"{'':>12} {'':>4} {'':>9} {'':>12} {'':>9} {'':>9} {'':>5} {'conc_win':>9} {'c_neg':>9} {'1-c_neg':>9} {'lift/(1-c)':>12}")
print("-"*130)

for (a, b, fourW) in [(14,44,2.0),(14,44,3.0),(14,44,4.0),(100,130,3.0),(100,130,4.0),(400,430,4.0)]:
    W = fourW/4
    result = analyze_negative_set(a, b, W)

    # Line 1: basic metrics
    print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {result['S_min']:>9.4f} {result['margin']:>+12.4e} {result['lift']:>9.4f} {result['neg_measure_frac']:>9.4f} {result['num_components']:>5} {result['S_neg_L1']:>9.4f}")
    # Line 2: concentration and correlation with lift
    one_minus_cneg = 1.0 - result['c_neg'] if result['c_neg'] < 1.0 else 0.0
    lift_over_one_minus_cneg = result['lift'] / one_minus_cneg if one_minus_cneg > 1e-10 else float('inf')
    print(f"{'':>12} {'':>4} {'':>9} {'':>12} {'':>9} {'':>9} {'':>5} {result['conc_window']:>9.4f} {result['c_neg']:>9.4f} {one_minus_cneg:>9.4f} {lift_over_one_minus_cneg:>12.4f}")

print("-"*130)
print("\nINTERPRETATION:")
print("  neg_frac    : fraction of window where S(r) < 0")
print("  #comp       : number of connected components of {S<0}")
print("  Sneg_L1     : total negative depth (integral of max(0, -S) over window)")
print("  conc_window : top eigenvalue of concentration onto full window (→1 = prolate saturates)")
print("  c_neg       : top eigenvalue of concentration onto {S<0} (can band-limited functions concentrate there?)")
print("  1-c_neg     : 'deficit' (if c_neg→1, functions concentrate fully; if c_neg<<1, they don't)")
print("  lift/(1-c)  : does the lift track the concentration deficit? If ~const or grows, uncertainty principle governs")
print("\nUNCERTAINTY HYPOTHESIS: If c_neg is bounded away from 1, band-limited functions of type 2·4W")
print("cannot fully concentrate on {S<0}. The margin compensates by the 'gap' (1−c_neg). This is why")
print("the lift (=margin−S_min) is non-trivial: the band-limiting constraint forces positivity.")
print("="*130)
