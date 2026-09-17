#!/usr/bin/env python3
"""
D-H certified-structure discriminator run (consumes the VALIDATED numeric archimedean block from
dh_arch_numeric.make_psiR). Builds:
   tau_DH  = -Q_arch(nu=3/4, cond=5)  -  Q_prime_DH(signed Lambda_D)      [NO W02: D-H is ENTIRE]
   tau_zeta= W02(pole) - Q_arch(nu=1/4,cond=1) - Q_prime_zeta(Lambda)      [Groskin, for side-by-side]
and reports min-eig / inertia for both on the SAME basis.
NOTE: the arch block is high-precision NUMERICAL (validated against Groskin's CERTIFIED closed form for zeta in
dh_arch_numeric.py) -> this run is high-precision numerical, NOT yet Arb-certified. Honest tier: strong numerical
evidence; Arb certification of the D-H arch is the follow-on. PREDICTION: min-eig(tau_DH) < 0 robustly (D-H has
off-line zeros) while min-eig(tau_zeta) ~ 0+ (razor).
"""
import sys, os, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "external", "groskin_2607.02828", "code"))
import mpmath as mp
mp.mp.dps = 25
import numpy as np
from dh_arch_numeric import make_psiR
from _dh_negative_control import lambda_D  # signed D-H generalized von Mangoldt


def vonmangoldt(N):
    lam = [0.0]*(N+1)
    for p in range(2, N+1):
        if lam[p] == 0:
            lp = math.log(p); pk = p
            while pk <= N:
                lam[pk] = lp; pk *= p
    return lam


def arch_matrix(N, Lval, nu, cond):
    psiR = make_psiR(Lval, nu, cond)
    vals = {k: psiR(k) for k in range(0, N+1)}
    h = mp.mpf(10)**(-6)
    diag = {n: (psiR(n+h)-psiR(n-h))/(2*h) for n in range(0, N+1)}
    A = np.zeros((2*N+1, 2*N+1))
    def val(k): return vals[abs(k)]
    for i in range(2*N+1):
        n = i-N
        for j in range(2*N+1):
            m = j-N
            if n == m:
                A[i, j] = float(diag[abs(n)])          # psi_R'(n) = -WR diagonal
            else:
                A[i, j] = float((val(m)-val(n))/(m-n))
    return A  # = -Q_arch (WR enters tau as -WR, and psiR already carries the -sign)


def prime_matrix(N, Lval, weights, positions):
    L = Lval
    A = np.zeros((2*N+1, 2*N+1))
    for i in range(2*N+1):
        n = i-N
        for j in range(2*N+1):
            m = j-N
            s = 0.0
            for w, y in zip(weights, positions):
                if n == m:
                    q = 2*(1-y/L)*math.cos(2*math.pi*n*y/L)
                else:
                    q = (math.sin(2*math.pi*m*y/L)-math.sin(2*math.pi*n*y/L))/(math.pi*(n-m))
                s += w*q
            A[i, j] = s
    return A


if __name__ == "__main__":
    c, N = 14, 10
    Lval = float(mp.log(c))
    # zeta prime side: Lambda(n) on prime powers
    VM = vonmangoldt(c)
    zpos = [math.log(n) for n in range(2, c+1) if VM[n] > 0]
    zw = [VM[n]/math.sqrt(n) for n in range(2, c+1) if VM[n] > 0]
    # D-H prime side: signed Lambda_D(n) on ALL n<=c
    LD = lambda_D(c)
    dpos = [math.log(n) for n in range(2, c+1)]
    dw = [LD[n]/math.sqrt(n) for n in range(2, c+1)]

    print(f"D-H vs zeta Weil min-eig on same basis (c={c}, N={N}); arch=validated-numeric, tier=high-prec numerical")
    archZ = arch_matrix(N, Lval, mp.mpf(1)/4, 1)   # -Q_arch zeta ; tau_zeta = W02 + archZ - primeZ ... see note
    archD = arch_matrix(N, Lval, mp.mpf(3)/4, 5)   # -Q_arch D-H (no pole)
    primeZ = prime_matrix(N, Lval, zw, zpos)
    primeD = prime_matrix(N, Lval, dw, dpos)

    # tau = W02(pole, zeta only) - WR - Wp ; archZ/archD already = -WR ; so tau = W02 + archX - primeX
    # zeta W02 (pole block) from Groskin closed form:
    from flint import arb, ctx
    ctx.prec = 200
    La = arb(c).log(); PI = arb.pi(); sp2 = 16*PI*PI; l2 = La*La
    pref02 = 32*La*(La/4).sinh()**2
    W02 = np.zeros((2*N+1, 2*N+1))
    for i in range(2*N+1):
        n = i-N
        for j in range(2*N+1):
            m = j-N
            num = l2 - sp2*m*n; den = (l2+sp2*m*m)*(l2+sp2*n*n)
            W02[i, j] = float((pref02*num/den).mid())

    tauZ = W02 + archZ - primeZ
    tauD = 0.0 + archD - primeD            # NO pole for entire D-H
    evZ = np.linalg.eigvalsh(0.5*(tauZ+tauZ.T))
    evD = np.linalg.eigvalsh(0.5*(tauD+tauD.T))
    print(f"  zeta : min-eig={evZ[0]:+.5e}  n_neg={int((evZ<-1e-9).sum())}  (expect ~0+ razor / PSD)")
    print(f"  D-H  : min-eig={evD[0]:+.5e}  n_neg={int((evD<-1e-9).sum())}  (expect <0: off-line zeros)")
    print(f"  => {'DISCRIMINATED: D-H min-eig robustly below zeta' if evD[0] < evZ[0]-1e-6 else 'inconclusive at this (c,N) — push scale'}")
