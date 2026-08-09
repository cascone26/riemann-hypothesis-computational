#!/usr/bin/env python3
"""§844: is the whole Weil form a MULTIPLICATION OPERATOR with symbol S(r)=(1/2π)Ω(r)−P(r), where
P(r)=(1/π)Σ_{log n<4W}Λ(n)/√n cos(log n·r) is the prime oscillation from the explicit formula? If so:
(a) margin = λmin(Q,G) ≥ min_window S(r)  (Rayleigh bound for a compressed multiplication operator),
(b) min_window S ≥ 0 ⟹ margin ≥ 0 — a PROVABLE sufficient condition (no RH input).
Test: compute S(r), its window-min, and compare to the measured margin at several (window,4W). Also check
whether margin ≈ min_window S (symbol-min) or margin ≫ min S (band-limiting lifts it, the hard regime)."""
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

def analyze(a, b, W, dr=0.004):
    fourW = 4*W; L = b-a
    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    # symbols
    arch_sym = OMEGA/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]
    P = np.zeros_like(RGRID)
    for u, w in zip(ln, lw):
        P += w*np.cos(u*RGRID)
    P *= (1.0/np.pi)                       # prime symbol (1/pi) sum Lambda/sqrt(n) cos(log n r)
    S = arch_sym - P                       # Weil symbol (ignoring pole, small)
    win = (RGRID >= a) & (RGRID <= b)
    S_min = float(S[win].min()); arch_min = float(arch_sym[win].min()); P_max = float(P[win].max())
    # operator margin (generalized min eig of Q=pole+arch-prime vs G)
    N = max(4, round(0.75*L*2*W/np.pi)); CEN = list(np.linspace(a, b, N))
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    Pm = np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G = (Pm @ Pm.T)*dr
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
    return S_min, arch_min, P_max, margin, N

print("="*104)
print("§844 Weil-symbol check: S(r)=(1/2π)Ω(r)−P(r); is margin ≥ min_window S? is min S≥0 the safe regime?")
print("="*104)
print(f"{'window':>14} {'4W':>5} {'archsym_min':>11} {'Psym_max':>9} {'S_min':>10} {'margin':>13} {'margin≥S_min?':>13} {'S_min≥0?':>9}")
for (a, b, fourW) in [(14,44,2.0),(14,44,3.0),(14,44,4.0),(100,130,3.0),(100,130,4.0),(400,430,4.0)]:
    W = fourW/4
    Smin, amin, Pmax, mg, N = analyze(a, b, W)
    print(f"{f'[{a},{b}]':>14} {fourW:>5.1f} {amin:>11.4f} {Pmax:>9.4f} {Smin:>10.4f} {mg:>+13.4e} {str(mg>=Smin-1e-6):>13} {str(Smin>=0):>9}")
print("="*104)
print("If margin ≥ S_min always AND (S_min≥0 ⟹ margin≥0): the Weil form is a compressed multiplication operator,")
print("min_window[(1/2π)Ω − P] ≥ 0 is a PROVABLE sufficient condition for positivity (both sides = explicit-formula symbols).")
print("Where S_min<0 but margin>0: band-limiting lifts it — the hard (RH) regime.")
print("="*104)
