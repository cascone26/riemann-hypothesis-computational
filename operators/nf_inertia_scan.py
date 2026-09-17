#!/usr/bin/env python3
"""§947 NUMBER-FIELD inertia of the actual arch+pole-prime Weil form Q (the ζ case).

Function-field §946 established: RH_ff <=> #neg(T_N)=0, a ROBUST integer (invariant to violation magnitude).
Now the number field. Under RH (true for all tested ζ zeros) the claim from §889 is: M_sigma=arch-prime has
exactly ONE negative eigenvalue at the corner, the rank-<=2 pole PSD-lifts it, so #neg(Q)=0. Here we:
  (1) VERIFY #neg(Q)=0 across many windows / bandwidths (the number-field index = 0, RH-consistent);
  (2) characterize the inertia of M_sigma alone (arch-prime, no pole) — how many negative dirs the pole must lift;
  (3) report the NULLITY of Q (its razor kernel dim) — the analog of the ff N-2g excess subspace.
If #neg(M_sigma) is a small robust integer (the pole rank bounds it) and #neg(Q)=0 everywhere, the number-field
index framing holds: RH <=> a signature/inertia condition, and the OPEN problem is a signature FORMULA (index
theorem) that gives 0 without diagonalizing — the Hodge-index analog function fields have and number fields lack.
HANDLE (own code). Reuses the verified mats() construction from prime_scaling_safety.py."""
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
NMAX = 5000; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

def mats(a, b, W, N, dr=0.004):
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
    def fc(c):
        y = W*(0.5j-c); y2 = W*(0.5j+c); return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in CEN]); pole = 2.0*np.real(np.outer(fp, fp))
    return pole+arch, arch, prime, G

def gen_inertia(M, G, tol=1e-9):
    """inertia of M w.r.t. metric G: solve G^{-1/2} M G^{-1/2}. Returns (#neg,#zero,#pos, eigs)."""
    gw, gv = np.linalg.eigh(G)
    keep = gw > gw.max()*1e-12
    B = gv[:, keep]/np.sqrt(gw[keep])
    ev = np.linalg.eigvalsh(B.T @ M @ B)
    scale = max(1.0, np.abs(ev).max())
    neg = int(np.sum(ev < -tol*scale)); pos = int(np.sum(ev > tol*scale))
    return neg, len(ev)-neg-pos, pos, ev

print("="*100)
print("§947 NUMBER-FIELD inertia of the Weil form Q=arch+pole-prime (ζ). Under RH expect #neg(Q)=0 everywhere.")
print("="*100)
print(f"{'window':>11} {'4W':>4} {'N':>3} | {'#neg(M_σ=arch-prime)':>20} {'#neg(Q, pole lifted)':>20} {'min eig Q':>13} {'nullity Q':>10}")
cfgs = [(14,44,2.0,0.7),(14,44,3.0,0.7),(14,44,4.0,0.7),(14,44,4.0,0.9),(14,44,5.0,0.7),
        (60,90,4.0,0.7),(120,150,4.0,0.7),(300,330,4.0,0.7),(14,34,4.0,0.9),(14,24,4.0,0.9)]
for (a,b,fourW,frac) in cfgs:
    W = fourW/4; L=b-a; Ncap = L*2*W/np.pi; N = max(4, round(frac*Ncap))
    AP, arch, prime, G = mats(a,b,W,N)
    Msig = arch - prime
    nneg_m, nz_m, np_m, _ = gen_inertia(Msig, G)
    nneg_q, nz_q, np_q, evq = gen_inertia(AP - prime, G)
    print(f"{f'[{a},{b}]':>11} {fourW:>4.1f} {N:>3} | {nneg_m:>20} {nneg_q:>20} {evq.min():>+13.3e} {nz_q:>10}")
print("="*100)
print("READ: #neg(Q)=0 across all windows => number-field index is 0 (RH-consistent). #neg(M_σ) is the number")
print("of negative directions the rank-<=2 pole must lift; if it stays small & bounded by pole rank, the index")
print("framing is exact. The OPEN problem becomes: a signature FORMULA computing #neg(Q)=0 (index theorem),")
print("the number-field analog of the Hodge-index theorem function fields have (§946) and number fields lack.")
print("="*100)
