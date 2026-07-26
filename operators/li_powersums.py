#!/usr/bin/env python3
"""§305: Li via zero power-sums. 1-(1-1/rho)^n = sum_{j=1}^n C(n,j)(-1)^{j-1} rho^{-j}, so
lambda_n = sum_{j=1}^n C(n,j)(-1)^{j-1} S_j,  S_j = sum_rho rho^{-j} (over nontrivial zeros).
Cross-checks §304's lambda_n a SECOND way, and exposes the S_j = the natural slot for arithmetic
(S_j = S_j^arch - S_j^prime, prime part from zeta'/zeta = von Mangoldt). Charisma: REPORT."""
import mpmath as mp
from mpmath import binomial as C
mp.mp.dps=30
gammas=[mp.mpf(x) for x in open("operators/_zeros150.txt").read().split()]
print(f"{len(gammas)} zeros",flush=True)
# S_j = sum over all nontrivial zeros rho^{-j}. Pair rho=1/2+i g and conj 1/2-i g.
def Sj(j):
    s=mp.mpf(0)
    for g in gammas:
        rho=mp.mpf(1)/2+1j*g
        s+=rho**(-j)+ (mp.conj(rho))**(-j)
    return s   # real (conjugate-symmetric)
J=8
S=[None]+[Sj(j) for j in range(1,J+1)]
print("zero power-sums S_j = sum_rho rho^-j (converge fast for j>=2):",flush=True)
for j in range(1,J+1):
    print(f"  S_{j} = {mp.nstr(mp.re(S[j]),10)}",flush=True)
# lambda_n via power sums (truncating j<=min(n,J)); compare to direct §304 values
print("lambda_n via S_j  vs  direct (§304, 150-zero partial):",flush=True)
direct={1:0.020642337,2:0.082532255,3:0.18555862,4:0.32953667,5:0.51420876,6:0.73924535}
for n in range(1,7):
    lam=sum(C(n,j)*(-1)**(j-1)*mp.re(S[j]) for j in range(1,n+1))
    print(f"  lambda_{n}: powersum={mp.nstr(lam,8)}  direct={direct[n]}  match={'YES' if abs(lam-direct[n])<1e-5 else 'no'}",flush=True)
# The arch/prime split lives in S_j. S_1 has a known closed value; higher S_j = arch consts - prime sums.
# Show the KNOWN arch constant for S_1 = 1 + gamma_E/2 - (1/2)ln(4pi)  (Riemann-von Mangoldt):
S1_arch=1+mp.euler/2-mp.log(4*mp.pi)/2
print(f"  S_1: computed={mp.nstr(mp.re(S[1]),8)}  known-closed 1+gammaE/2-ln(4pi)/2={mp.nstr(S1_arch,8)}  match={'YES' if abs(mp.re(S[1])-S1_arch)<1e-3 else 'partial(tail)'}",flush=True)
print("DONE",flush=True)
