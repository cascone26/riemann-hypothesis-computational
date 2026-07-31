#!/usr/bin/env python3
"""Baez-Duarte distance d_N for the Nyman-Beurling system (q=1) — likely-unpublished numerics.
d_N^2 = C - F^T G^{-1} F  (optimal LS),  C = 1,  G[u,v]=A(u,v)/u (validated §477),
F_n = (1/2pi i) int_{Re s=1/2} n^{-s} zeta(s)/(s(1-s)) ds = n^{-1/2}(1/2pi) int n^{-it} zeta(1/2+it)/(1/4+t^2) dt.
Also the Baez-Duarte Mobius coefficients lambda_{n,N}=mu(n)(1-log n/log N): d_N^2(lambda)=C-2 sum lam F + lam^T G lam.
Conjecture (Baez-Duarte-Balazard-Landreau-Saias): d_N^2 ~ C_BD/log N, C_BD = 2 + gamma - log(4 pi) ~ 0.0463.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from nb_gram import build_gram
mp.mp.dps = 20

C_BD = float(2 + mp.euler - mp.log(4*mp.pi))

def F_n(n, T=600.0):
    ln = mp.log(n)
    def integrand(t):
        z = mp.zeta(mp.mpf('0.5')+1j*t)
        return (mp.e**(-1j*ln*t) * z / (mp.mpf('0.25')+t*t)).real
    I = mp.quad(integrand, [-T,-100,-30,-10,0,10,30,100,T])
    return float(mp.mpf(n)**mp.mpf('-0.5') * I/(2*mp.pi))

def mobius(N):
    mu = np.ones(N+1); is_c = np.zeros(N+1, bool); primes=[]
    mu[0]=0
    small = np.zeros(N+1, int)
    # simple mobius sieve
    m = np.ones(N+1, dtype=int); m[0]=0
    prime = np.ones(N+1, bool); prime[:2]=False
    mo = np.ones(N+1, dtype=int)
    for p in range(2, N+1):
        if prime[p]:
            for k in range(p, N+1, p): prime[k]=False if k>p else prime[k]
    # redo properly
    mo = np.ones(N+1, dtype=int); mo[0]=0
    least = np.zeros(N+1, dtype=int)
    for i in range(2, N+1):
        if least[i]==0:
            for k in range(i, N+1, i):
                if least[k]==0: least[k]=i
    for i in range(2, N+1):
        x=i; last=0; sqfree=True; cnt=0
        while x>1:
            p=least[x]; e=0
            while x%p==0: x//=p; e+=1
            if e>1: sqfree=False; break
            cnt+=1
        mo[i] = 0 if not sqfree else (-1)**cnt
    mo[1]=1
    return mo

if __name__ == "__main__":
    Nmax = int(sys.argv[1]) if len(sys.argv)>1 else 60
    print(f"Computing F_n for n=1..{Nmax} (mpmath integrals)...", flush=True)
    F = np.array([F_n(n) for n in range(1, Nmax+1)])
    print(f"F_1..5 = {np.round(F[:5],5)}")
    mo = mobius(Nmax)
    print(f"\n{'N':>4} {'d_N^2 (optimal)':>16} {'d_N^2 (Mobius)':>15} {'C_BD/logN':>11} {'opt*logN':>10}")
    for N in [5,10,20,30,40,50,60,80,100][:99]:
        if N>Nmax: break
        M = build_gram(N); FN = F[:N]
        x = np.linalg.solve(M, FN)
        d2_opt = 1.0 - FN@x
        lam = mo[1:N+1]*(1 - np.log(np.arange(1,N+1))/np.log(N))
        d2_mob = 1.0 - 2*lam@FN + lam@M@lam
        print(f"{N:>4} {d2_opt:>16.6e} {d2_mob:>15.6e} {C_BD/np.log(N):>11.4e} {d2_opt*np.log(N):>10.4f}", flush=True)
    print(f"\nC_BD = 2+gamma-log(4pi) = {C_BD:.6f}. If RH & conjecture: d_N^2*logN -> C_BD.")
    print("opt <= mob (optimal beats fixed Mobius coeffs). d_N^2 decreasing => consistent with RH.")
