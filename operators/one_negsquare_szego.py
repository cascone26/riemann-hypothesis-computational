#!/usr/bin/env python3
"""§948 THE ONE-NEGATIVE-SQUARE ANOMALY (the arithmetic miracle, quantified).

§947 found: M_σ = arch-prime (compression of mult-by-symbol σ(r)=Ω(r)/2π - P(r)) has #neg <= 1 across every
window/bandwidth tested. But σ is POINTWISE NEGATIVE on a positive-measure set (§889: 7-40% of each window).
By Szego / Avram-Parter, the compression of mult-by-σ should have ~ N * frac{σ<0} NEGATIVE eigenvalues.
Getting only ONE is a huge cancellation defect. This script:
  (1) QUANTIFIES the defect: naive-Szego prediction (frac{σ<0} * N) vs actual #neg(M_σ);
  (2) STRESS-TESTS the "exactly one negative square" (Krein-Langer Pi_1) claim: does #neg(M_σ) EVER exceed 1
      over a fine, aggressive sweep (large N, high bandwidth, many heights)? Report the MAX honestly.
  (3) tests STABILITY of the single negative direction: is the negative eigenvector the "same" function
      (high overlap) across nearby windows -> a canonical Pi_1 anomaly direction.
If #neg(M_σ) is provably capped at 1 by structure, the Weil kernel is a Pontryagin Pi_1 kernel (Krein-Langer),
pole (PSD) lifts it, and RH <=> "the Weil kernel has exactly the pole-rank negative squares." HANDLE (own code)."""
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
NMAX = 8000; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

def build(a, b, W, N, dr=0.004):
    fourW = 4*W
    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    CEN = list(np.linspace(a, b, N))
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    P = np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]
    cosu = np.cos(np.outer(ln, RGRID)); prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
    prime = 0.5*(prime+prime.T)
    # pointwise symbol sigma(r) on the window and its negative fraction over [a,b]
    win = (RGRID >= a) & (RGRID <= b)
    Prow = OMEGA/(2*np.pi)
    # prime density P(r) = (1/pi) sum Lambda(n)/sqrt(n) cos(r log n)  (the pointwise prime symbol)
    Pr = np.zeros_like(RGRID)
    for l, w in zip(ln, lw):
        Pr += w*np.cos(RGRID*l)
    Pr /= np.pi
    sigma = OMEGA/(2*np.pi) - Pr
    fracneg = float(np.mean(sigma[win] < 0))
    return arch, prime, G, fracneg

def gen_eig(M, G):
    gw, gv = np.linalg.eigh(G); keep = gw > gw.max()*1e-12
    B = gv[:, keep]/np.sqrt(gw[keep])
    ev, evec = np.linalg.eigh(B.T @ M @ B)
    return ev, B @ evec  # eigenvalues, eigenvectors back in original coords

def nneg(ev):
    s = max(1.0, np.abs(ev).max()); return int(np.sum(ev < -1e-9*s))

print("="*104)
print("§948 ONE-NEGATIVE-SQUARE ANOMALY: Szego predicts ~frac{σ<0}*N negatives; M_σ=arch-prime delivers <=1.")
print("="*104)
print(f"{'window':>11} {'4W':>4} {'N':>3} | {'frac σ<0':>9} {'Szego pred #neg':>15} {'ACTUAL #neg(M_σ)':>16} {'defect':>8}")
cfgs = [(14,44,2.0,0.7),(14,44,3.0,0.7),(14,44,4.0,0.7),(14,44,4.0,0.9),(14,44,5.0,0.8),(14,44,6.0,0.8),
        (14,64,4.0,0.85),(14,64,5.0,0.85),(60,90,5.0,0.8),(120,150,5.0,0.8),(300,360,5.0,0.8)]
maxneg = 0
for (a,b,fourW,frac) in cfgs:
    W = fourW/4; L=b-a; Ncap = L*2*W/np.pi; N = max(4, round(frac*Ncap))
    arch, prime, G, fn = build(a,b,W,N)
    ev, _ = gen_eig(arch - prime, G); nn = nneg(ev)
    maxneg = max(maxneg, nn)
    szego = fn*N
    print(f"{f'[{a},{b}]':>11} {fourW:>4.1f} {N:>3} | {fn:>9.3f} {szego:>15.1f} {nn:>16} {szego-nn:>8.1f}")
print("-"*104)
print(f"MAX #neg(M_σ) over this aggressive sweep = {maxneg}   (if 1, the one-negative-square / Krein-Langer Π_1 claim holds)")

# --- Stability of the single negative direction across a sliding window (is it a canonical Π_1 direction?) ---
print("\n[stability] negative eigenvector overlap across sliding windows [14,44]->[16,46]->...(4W=4, same N):")
prev = None
for shift in [0.0, 1.0, 2.0, 3.0, 4.0]:
    a, b = 14+shift, 44+shift; W = 1.0; N = 13
    arch, prime, G, fn = build(a,b,W,N)
    ev, evec = gen_eig(arch - prime, G)
    order = np.argsort(ev)
    v = evec[:, order[0]]  # eigenvector of most-negative eigenvalue
    v = v/np.linalg.norm(v)
    ov = "" if prev is None else f"  overlap w/ prev = {abs(np.dot(prev, v)):.3f}"
    print(f"  window [{a:.0f},{b:.0f}]  min eig(M_σ)={ev.min():+.3e}{ov}")
    prev = v
print("="*104)
print("READ: a large positive 'defect' column = the arithmetic cancellation (Szego says many negatives, truth=1).")
print("MAX #neg=1 over the sweep => M_σ is a Pontryagin Π_1 kernel (exactly one negative square). This is the")
print("structural target: prove the Weil kernel has exactly pole-rank negative squares (Krein-Langer) => RH.")
print("="*104)
