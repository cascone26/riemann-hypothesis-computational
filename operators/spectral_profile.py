#!/usr/bin/env python3
"""§850: is the Weil margin an ISOLATED small eigenvalue (positivity governed by ONE mode → rank-one
structure, connecting to CCM's rank-one-perturbation operator) or part of a cluster (genuinely
high-dimensional)? Compute the FULL generalized spectrum of Q=pole+arch−prime vs G at the corner and a
few windows; report λ1 (margin), λ2, the gap λ2/λ1, and how many eigenvalues are 'small' (< 0.05).
If λ2/λ1 is huge (margin isolated) → the obstruction is effectively rank-one; if λ's cluster near 0 →
high-dimensional. Float (diagnostic)."""
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

def spectrum(a, b, W, dr=0.004):
    fourW = 4*W; L = b-a
    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    N = max(4, round(0.75*L*2*W/np.pi)); CEN = list(np.linspace(a, b, N))
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    Pm = np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G = (Pm @ Pm.T)*dr
    arch = (Pm*OMEGA[None, :]) @ Pm.T * dr/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]
    cosu = np.cos(np.outer(ln, RGRID)); prime = np.zeros((N, N))
    for j in range(N):
        gp = ((Pm[j][None, :]*Pm) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
    prime = 0.5*(prime+prime.T)
    def fc(c):
        y = W*(0.5j-c); y2 = W*(0.5j+c); return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in CEN]); pole = 2.0*np.real(np.outer(fp, fp))
    w = eigh(pole+arch-prime, G, eigvals_only=True)   # ascending generalized eigenvalues
    return w, N

print("="*96)
print("§850 Weil form spectral profile — is the margin ISOLATED (rank-one) or CLUSTERED?")
print("="*96)
print(f"{'window':>12} {'4W':>4} {'N':>3} {'λ1(margin)':>12} {'λ2':>10} {'λ2/λ1':>10} {'#small(<0.05)':>13} {'λ_top':>8}")
for (a, b, fourW) in [(14,44,3.0),(14,44,4.0),(14,44,2.0),(100,130,4.0),(400,430,4.0)]:
    W = fourW/4; w, N = spectrum(a, b, W)
    l1, l2 = w[0], w[1]
    nsmall = int((w < 0.05).sum())
    ratio = l2/l1 if l1 > 0 else float('inf')
    print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {N:>3} {l1:>+12.4e} {l2:>+10.4e} {ratio:>10.2e} {nsmall:>13} {w[-1]:>8.3f}")
print("="*96)
print("λ2/λ1 huge + #small=1 → margin isolated, positivity is effectively RANK-ONE (→ CCM rank-one operator).")
print("λ's cluster near 0 (#small large) → high-dimensional obstruction (harder).")
print("="*96)
