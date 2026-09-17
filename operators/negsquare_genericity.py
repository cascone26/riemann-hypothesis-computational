#!/usr/bin/env python3
"""§949 IS THE ONE-NEGATIVE-SQUARE PROPERTY ARITHMETIC OR STRUCTURAL? (the make-or-break disentangle)

§948 found M_σ = arch-prime has exactly ONE negative square (Π_1) for the REAL primes. Decisive question:
does κ=1 come from the ARITHMETIC (real von Mangoldt / real log-primes — genuine RH-content), or is it a
STRUCTURAL fact of band-limited compressions of any such symbol (generic — would mean 'prove Π_1' needs no
arithmetic, and the RH-content must live entirely in the pole-lift)?

Build M = arch - prime_variant and count negative squares κ = #neg(M) for:
  (0) REAL: true von Mangoldt weights Λ(n)/√n at true log n            [baseline, expect κ=1]
  (P) PERTURBED log-primes: log n -> log n * (1+δ ξ_n), ξ_n ~ U(-1,1)  [§944-style arithmetic break]
  (S) SCRAMBLED frequencies: random freqs in same range, real weights   [destroys prime arithmetic entirely]
  (W) PERMUTED weights: real freqs, weights shuffled across n            [right frequencies, wrong amplitudes]
  (R) RANDOM symbol: replace prime density by a random smooth field of matched negativity fraction
For each: report κ over increasing bandwidth. INTERPRETATION printed at end, honestly, both ways.
HANDLE (own code). arch is fixed (real archimedean); only the prime/symbol side varies."""
import numpy as np
from scipy.special import digamma

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam
NMAX = 8000; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

def build(a, b, W, N, freqs, weights, dr=0.004):
    """M = arch - prime_variant vs metric G, using given (freqs, weights) for the prime side."""
    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    CEN = list(np.linspace(a, b, N))
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    P = np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    cosu = np.cos(np.outer(freqs, RGRID)); prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*weights[None, :]).sum(axis=1)
    prime = 0.5*(prime+prime.T)
    return arch, prime, G

def kappa(arch, prime, G):
    gw, gv = np.linalg.eigh(G); keep = gw > gw.max()*1e-12
    B = gv[:, keep]/np.sqrt(gw[keep]); ev = np.linalg.eigvalsh(B.T @ (arch-prime) @ B)
    s = max(1.0, np.abs(ev).max()); return int(np.sum(ev < -1e-9*s))

rng = np.random.default_rng(2026)
a, b = 14, 44
print("="*100)
print("§949 negative-square count κ=#neg(arch-prime) — REAL vs perturbed/scrambled arithmetic. Window [14,44].")
print("="*100)
print(f"{'4W':>4} {'N':>3} | {'(0)REAL':>8} {'(P)δ=1e-3':>10} {'(P)δ=1e-2':>10} {'(S)scramb':>10} {'(W)permwt':>10} {'(R)rand':>9}")
for fourW in [3.0, 4.0, 5.0, 6.0]:
    W = fourW/4; L=b-a; Ncap = L*2*W/np.pi; N = max(4, round(0.8*Ncap))
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]; m = len(ln)
    # (0) real
    A,Pm,G = build(a,b,W,N,ln,lw); k0 = kappa(A,Pm,G)
    # (P) perturbed log-primes, two deltas
    kP=[]
    for delta in [1e-3, 1e-2]:
        lnp = ln*(1+delta*rng.uniform(-1,1,m)); A,Pm,G = build(a,b,W,N,lnp,lw); kP.append(kappa(A,Pm,G))
    # (S) scrambled: random freqs in [min ln, max ln], real weights
    fs = rng.uniform(ln.min(), ln.max(), m); A,Pm,G = build(a,b,W,N,fs,lw); kS = kappa(A,Pm,G)
    # (W) permuted weights at real freqs
    wp = rng.permutation(lw); A,Pm,G = build(a,b,W,N,ln,wp); kW = kappa(A,Pm,G)
    # (R) random symbol: random freqs + random positive weights of matched total mass
    fr = rng.uniform(0.1, ln.max(), m); wr = rng.uniform(0,1,m); wr *= lw.sum()/wr.sum()
    A,Pm,G = build(a,b,W,N,fr,wr); kR = kappa(A,Pm,G)
    print(f"{fourW:>4.1f} {N:>3} | {k0:>8} {kP[0]:>10} {kP[1]:>10} {kS:>10} {kW:>10} {kR:>9}")
print("="*100)
print("INTERPRET (both ways, honestly):")
print(" • If REAL gives κ=1 but perturbed/scrambled give κ>1  => Π_1 is ARITHMETIC: real primes are special,")
print("   'prove κ=1' is genuine RH-content (the target lane is real and non-circular-ish).")
print(" • If perturbed/scrambled ALSO give κ=1 (or ≤1)         => Π_1 is STRUCTURAL/generic: band-limiting caps")
print("   negative squares regardless of arithmetic; then κ=1 alone can't imply RH — the RH-content lives ENTIRELY")
print("   in the POLE-LIFT (does the rank-≤2 pole lift THIS particular negative direction to ≥0), not in the count.")
print("="*100)
