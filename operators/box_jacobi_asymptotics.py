#!/usr/bin/env python3
"""
cycle 332 — asymptotic LAW of the box-Jacobi (candidate heig) coefficients from the zero density N(T).

box measure dnu = sum_gamma x_gamma delta_{x_gamma}, x_gamma = 1/(1/4+gamma^2), weight x_gamma.
Its Jacobi matrix J_nu (Favard) is the candidate heig (J_nu^{-1} spectrum = {1/4+gamma^2} = boxEigen).
The INTEGRATED DENSITY OF STATES is UNCONDITIONAL:
    n(E) = #{gamma : x_gamma >= E} = #{gamma : gamma <= T(E)} = N(T(E)),  T(E)=sqrt(1/E - 1/4),
where N(T) = (T/2pi) log(T/2pi) - T/2pi + 7/8 + ...  (Riemann-von Mangoldt, KNOWN unconditionally).

SEMICLASSICAL/WEYL PREDICTION for the Jacobi diagonal: a_n ~ E_n, where N(T(E_n)) = n, i.e.
    T_n = N^{-1}(n),   a_n^pred = 1/(1/4 + T_n^2).
This is a law expressed PURELY through N(T) -- no zero locations, unconditional.

TEST: compute the TRUE (a_n,b_n) via Lanczos on 2000 zeros (stable), compare a_n to the N(T)-prediction,
and fit the decay of b_n. If a_n ~ 1/(1/4+T_n^2) cleanly, the operator's structure is N(T)-pinned.
Verify-before-ship: this uses zeros to DISCOVER/verify the law; the LAW itself (via N(T)) is unconditional.
"""
import numpy as np
from mpmath import mp, mpf, findroot, log as mlog, pi as mpi

mp.dps = 20

# --- true Jacobi coeffs via Lanczos on zeros ---
G = np.load('/Users/scones/projects/riemann/operators/zeros2000.npy').astype(float)
x = 1.0/(0.25 + G**2)          # nodes
w = x/ x.sum()                 # dnu = x dmu ; normalized
def lanczos(nodes, wt, K):
    D = len(nodes); V = np.zeros((D, K)); a = np.zeros(K); b = np.zeros(K)
    q = np.sqrt(wt); q = q/np.linalg.norm(q); qp = np.zeros(D); bprev = 0.0
    for k in range(K):
        V[:, k] = q
        v = nodes*q; aa = float(q@v); a[k] = aa
        v = v - aa*q - bprev*qp
        v = v - V[:, :k+1] @ (V[:, :k+1].T @ v)   # FULL reorth
        v = v - V[:, :k+1] @ (V[:, :k+1].T @ v)   # twice (kills ghosts)
        nb = float(np.linalg.norm(v))
        if nb < 1e-14:
            a = a[:k+1]; b = b[:k+1]; break
        if k+1 < K: b[k+1] = nb
        qp = q; q = v/nb; bprev = nb
    return a, b
K = 50
a, b = lanczos(x, w, K)

# --- N(T) Riemann-von Mangoldt and its inverse ---
TWO_PI = 2*np.pi
def Nsmooth(T):
    # real float; monotone increasing for T > 2pi
    if T <= TWO_PI: return -0.125 - (TWO_PI - T)  # dummy decreasing continuation (keeps bracket sign)
    return T/TWO_PI*np.log(T/TWO_PI) - T/TWO_PI + 7.0/8.0
def Ninv(nn):
    # bracketed bisection for Nsmooth(T)=nn, T in (2pi, hi]; valid for nn>=~1
    lo, hi = TWO_PI + 1e-9, 20.0
    while Nsmooth(hi) < nn: hi *= 1.5
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if Nsmooth(mid) < nn: lo = mid
        else: hi = mid
    return 0.5*(lo+hi)

# --- compare a_n to the unconditional N(T)-prediction a_n^pred = 1/(1/4+T_n^2), T_n=Ninv(n+1/2) ---
print(f"{'n':>3} {'a_n(true)':>13} {'a_n^pred(N(T))':>15} {'ratio':>8} {'b_n(true)':>13}")
Nfit = 40
avals=[]; pvals=[]
for n in range(Nfit):
    Tn = Ninv(n + 0.5)               # midpoint convention: n+1/2 zeros below
    apred = 1.0/(0.25 + Tn**2)
    ratio = a[n]/apred if apred>0 else float('nan')
    avals.append(a[n]); pvals.append(apred)
    if n < 20 or n % 5 == 0:
        print(f"{n:>3} {a[n]:>13.6e} {apred:>15.6e} {ratio:>8.4f} {b[n]:>13.6e}")

# --- fit power laws (log-log slope) for a_n, b_n over a mid-range window (avoid n=0 transient & trunc tail) ---
ns = np.arange(5, 35)
la = np.log(np.array([a[n] for n in ns])); lb = np.log(np.array([b[n] for n in ns])); ln = np.log(ns.astype(float))
sa = np.polyfit(ln, la, 1)[0]; sb = np.polyfit(ln, lb, 1)[0]
print(f"\npower-law slopes over n in [5,35): a_n ~ n^{sa:.3f} ,  b_n ~ n^{sb:.3f}")
print(f"  (Weyl heuristic: T_n~2pi n/log n => 1/(1/4+T_n^2) ~ (log n)^2/n^2 => slope ~ -2 with log correction)")

# --- ratio a_n/a_n^pred: is it converging to a constant (=> law holds) or drifting? ---
r = np.array([avals[n]/pvals[n] for n in range(Nfit)])
print(f"\na_n/a_n^pred: n=10..39 mean={r[10:].mean():.4f} std={r[10:].std():.4f}  (const => N(T)-law holds)")
print(f"  tail ratios n=30..39: {np.round(r[30:40],4)}")
