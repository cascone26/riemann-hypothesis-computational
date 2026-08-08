#!/usr/bin/env python3
"""§832-VERIFY (float, fast) — is the 4W=2.2 knife-edge margin collapse REAL prime-side thinning or a
bad-span conditioning artifact? Float margin (generalized eig of Q,G) + λmin(G) + cond(G), sweeping
N at each 4W. Decisive: at 4W=2.2 does the margin stay small even at conservative WELL-conditioned N
(→ real thinning) or recover to ~0.1 when G is healthy (→ conditioning)? The arb-CERTIFIED margin
values already exist (§832 ladder: 0.151/0.145/0.103/5.6e-5); this only diagnoses the mechanism.
Faithful float build = operators/l2_fejer.build() with W parametrized."""
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
NMAX = 300; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)
RGRID = np.arange(-200.0, 200.0, 0.0025); dr = 0.0025
OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)

def fej(x, W):
    y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
    out[nz] = (np.sin(y[nz])/y[nz])**2; return out
def Frow(c, r, W): return fej(r-c, W) + fej(r+c, W)

def build(cs, W):
    n = len(cs); fourW = 4*W
    P = np.array([Frow(c, RGRID, W) for c in cs])
    G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]
    prime = np.zeros((n, n))
    if ln.size:
        cosu = np.cos(np.outer(ln, RGRID))
        for j in range(n):
            gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
            prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
        prime = 0.5*(prime+prime.T)
    def fc(c):
        y = W*(0.5j-c); y2 = W*(0.5j+c)
        return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in cs]); pole = 2.0*np.real(np.outer(fp, fp))
    return pole+arch-prime, G

def diag(cs, W):
    Q, G = build(cs, W)
    gw = np.linalg.eigvalsh(G); lg = float(gw[0]); cg = float(gw[-1]/max(gw[0], 1e-300))
    m = float(eigh(Q, G, eigvals_only=True)[0])
    return m, lg, cg

print("="*100)
print("§832-VERIFY (float) — knife-edge margin: real prime thinning vs conditioning?")
print("="*100)
print(f"{'4W':>6} {'Ncap':>6} {'N':>4} {'N/Ncap':>7} {'nprime':>7} {'float margin':>14} {'λmin(G)':>10} {'cond(G)':>10}  health")
for fourW in [1.60, 1.90, 2.20, 2.50]:
    W = fourW/4; Ncap = 15*fourW/np.pi; npr = int((LOGN < fourW).sum())
    Ns = sorted(set(max(4, round(f*Ncap)) for f in [0.45, 0.55, 0.65, 0.76]))
    print("-"*100)
    for N in Ns:
        CEN = list(np.linspace(14, 44, N)); m, lg, cg = diag(CEN, W)
        health = "G-healthy" if (lg > 0.5 and cg < 30) else ("G-loading" if cg < 1e3 else "G-DEGEN")
        print(f"{fourW:>6.2f} {Ncap:>6.2f} {N:>4} {N/Ncap:>7.2f} {npr:>7} {m:>+14.5e} {lg:>10.3f} {cg:>10.1f}  {health}")
print("="*100)
print("READ: 4W=2.2/2.5 margin small even at conservative N with G-healthy → REAL prime-side thinning (knife-edge).")
print("      margin recovers to ~0.1 at well-conditioned N → the sharp collapse was conditioning.")
print("="*100)
