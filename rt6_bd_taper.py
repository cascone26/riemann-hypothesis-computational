#!/usr/bin/env python3
"""§437 (Báez-Duarte mechanism): the near-optimal c_n ≈ μ(n)·τ(n/N) — what is the TAPER τ? The damping profile
is the actual mechanism making d_N→0 (Σμ(n)γ_n=1 diverges; the taper regularizes it). Compute c_n/μ(n) for
squarefree n, plot vs n/N. Conjectured shape ~ (1 − log n/log N) (linear in log). Characterize it."""
import numpy as np
from sympy import mobius
J=400000; j=np.arange(1,J+1,dtype=np.float64); w=1.0/(j*(j+1))
N=120
V=np.array([(j%n)/n for n in range(2,N+1)]); Vw=V*w[None,:]
G=Vw@V.T; b=-(Vw@np.ones(J)); c=np.linalg.solve(G,b)
ns=np.arange(2,N+1); mu=np.array([int(mobius(n)) for n in ns])
sf=mu!=0                       # squarefree
ratio=np.where(sf, c/np.where(mu!=0,mu,1), np.nan)
print(f"§437 Báez-Duarte taper c_n/μ(n) vs n/N (N={N}). squarefree n only.")
print(f"{'n':>4} {'n/N':>5} {'c_n/μ(n)':>9} {'1-logn/logN':>12}")
for i in range(0,len(ns),8):
    if sf[i]:
        print(f"{ns[i]:>4} {ns[i]/N:5.2f} {ratio[i]:9.4f} {1-np.log(ns[i])/np.log(N):12.4f}")
# fit c_n/μ(n) vs 1 - log n / log N
x=1-np.log(ns[sf])/np.log(N); y=ratio[sf]
A=np.polyfit(x,y,1)
r=np.corrcoef(x,y)[0,1]
print(f"\nfit c_n/μ(n) = {A[0]:.3f}·(1−logn/logN) + {A[1]:.3f}   (corr={r:.3f})")
# also vs plain n/N
x2=ns[sf]/N; A2=np.polyfit(x2,y,1); r2=np.corrcoef(x2,y)[0,1]
print(f"fit vs (1−n/N): corr with n/N = {np.corrcoef(1-x2,y)[0,1]:.3f}")
print("""
reading: if c_n/μ(n) ≈ (1 − log n/log N) with high corr, the taper is LOGARITHMIC — the Báez-Duarte mechanism
damps μ(n) linearly in log n up to the cutoff (the 'smoothed Möbius' that makes Σμγ=1 converge). This is the
concrete shape of the approximation-to-RH mechanism, computed. High corr = clean law.
""")
