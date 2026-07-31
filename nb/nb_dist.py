#!/usr/bin/env python3
"""Baez-Duarte distance d_N, ALL closed-form (validated §477/§477c):
  d_N^2 = C - F^T G^{-1} F,  C = 1,  G[u,v] = A(u,v)/u,  F_n = (gamma - 1 - log n)/n.
No zeta integrals. Conjecture (BDBLS): d_N^2 ~ C_BD/log N, C_BD = 2+gamma-log(4pi) ~ 0.0463.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nb_gram import build_gram
GAMMA = 0.5772156649015329
C_BD = 2 + GAMMA - np.log(4*np.pi)

def F_vec(N):
    n = np.arange(1, N+1)
    return (GAMMA - 1 - np.log(n))/n

def mobius(N):
    least = np.zeros(N+1, dtype=int)
    for i in range(2, N+1):
        if least[i] == 0:
            for k in range(i, N+1, i):
                if least[k] == 0: least[k] = i
    mo = np.ones(N+1, dtype=int); mo[0] = 0
    for i in range(2, N+1):
        x = i; cnt = 0; sq = True
        while x > 1:
            p = least[x]; e = 0
            while x % p == 0: x //= p; e += 1
            if e > 1: sq = False; break
            cnt += 1
        mo[i] = 0 if not sq else (-1)**cnt
    return mo

if __name__ == "__main__":
    Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    Gfull = build_gram(Nmax); Ffull = F_vec(Nmax); mo = mobius(Nmax)
    print(f"Baez-Duarte d_N (all closed-form). C_BD = 2+gamma-log(4pi) = {C_BD:.6f}")
    print(f"{'N':>5} {'d_N^2 (optimal)':>16} {'d_N^2 (Mobius)':>15} {'opt*logN':>10} {'mob*logN':>10} {'cond(G)':>10}")
    for N in [5,10,20,40,80,120,160,240,320,400]:
        if N > Nmax: break
        G = Gfull[:N,:N]; F = Ffull[:N]
        x = np.linalg.solve(G, F)
        d2_opt = 1.0 - F@x
        lam = mo[1:N+1]*(1 - np.log(np.arange(1,N+1))/np.log(N))
        d2_mob = 1.0 - 2*lam@F + lam@G@lam
        cond = np.linalg.cond(G)
        print(f"{N:>5} {d2_opt:>16.6e} {d2_mob:>15.6e} {d2_opt*np.log(N):>10.4f} {d2_mob*np.log(N):>10.4f} {cond:>10.2e}", flush=True)
    print(f"\nIf RH & BDBLS conjecture: d_N^2*logN -> C_BD = {C_BD:.5f}. Optimal <= Mobius. Monotone decrease => RH-consistent.")
