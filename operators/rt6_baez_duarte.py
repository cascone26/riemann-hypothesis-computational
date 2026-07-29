#!/usr/bin/env python3
"""§433 (NEW RH-equivalent, orthogonal to Weil): Nyman-Beurling-Báez-Duarte. RH ⟺ d_N→0 where
d_N²=min_c ||1 − Σ_{n≤N} c_n γ_n||², γ_n(t)=⌊t/n⌋−⌊t⌋/n ∈ L²([1,∞),dt/t²). KEY: γ_n is a STEP function
= −{j/n} on [j,j+1), so EXACTLY: G_{n,m}=Σ_j {j/n}{j/m}/(j(j+1)), b_n=−Σ_j {j/n}/(j(j+1)), ⟨1|1⟩=1.
d_N²=1−bᵀG⁻¹b; optimal c=G⁻¹b. Compute d_N vs N (does it →0? rate ~1/log N under RH?), and compare optimal
c_n to μ(n) (Báez-Duarte: Σμ(n)γ_n=1). Genuinely different lens; may reveal structure Weil positivity hides."""
import numpy as np
from sympy import mobius
J=300000
j=np.arange(1,J+1)
w=1.0/(j*(j+1))                       # weight; Σw = 1 (telescoping)
print(f"§433 Báez-Duarte. J={J} (tail~{1/(4*J):.1e}). check Σw={w.sum():.6f} (→1)")
def frac(n): return (j % n)/n        # {j/n} exactly for integer n
# build V rows for n=2..Nmax once
Nmax=70
V=np.array([frac(n) for n in range(2,Nmax+1)])   # (Nmax-1) x J
Vw=V*w[None,:]
Gfull=Vw@V.T                          # full Gram for n=2..Nmax
bfull=-(Vw@np.ones(J))                # b_n
onenorm=w.sum()                       # <1|1>
print(f"{'N':>4} {'d_N^2':>12} {'d_N':>10} {'d_N^2*logN':>11}  top optimal c_n (n=2..7)")
mu=np.array([int(mobius(n)) for n in range(2,Nmax+1)])
for N in [5,10,15,20,30,40,50,60,70]:
    k=N-1
    G=Gfull[:k,:k]; b=bfull[:k]
    c=np.linalg.solve(G,b)
    d2=onenorm - b@c
    d2=max(d2,0.0)
    print(f"{N:>4} {d2:12.6e} {np.sqrt(d2):10.5f} {d2*np.log(N):11.5f}  "
          f"{np.array2string(c[:6],precision=3,floatmode='fixed')}",flush=True)
# compare optimal c to mu(n) for a fixed N
N=40; k=N-1
c=np.linalg.solve(Gfull[:k,:k],bfull[:k])
print(f"\noptimal c_n vs μ(n), N={N}:")
for i,n in enumerate(range(2,12)):
    print(f"  n={n}: c_n={c[i]:+.4f}  μ(n)={mu[i]:+d}  c_n·? (BÁez-Duarte predicts c_n≈μ(n) near-optimal)")
print("""
reading: d_N² DECREASING toward 0 => consistent with RH (can't prove the limit). d_N²·log N →const => the
conjectured ~1/log N rate. optimal c_n tracking μ(n) => the Möbius structure (Σμ(n)γ_n=1) is the mechanism.
This is a fully independent RH-equivalent computed EXACTLY (step-function γ_n) — a clean cross-check lens.
""")
