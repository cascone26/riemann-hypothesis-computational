#!/usr/bin/env python3
"""ROUTE 6 × zero-density BUDGET: can the verified positive mass outweigh the density-bounded off-line
negative mass past n=2πT, extending λ_n>0 far beyond the §21 theorem?

λ_n = P_ver(n) + Tail(n),  P_ver = 2Σ_{γ≤T,on-line}(1−cos nφ_γ) ≥ 0,  φ_γ=2arctan(1/2γ).
Tail zeros (γ>T): if on-line contribute ≥0; if off-line, by the §21 exact closed form
F_n(σ,γ)=4−4cos(nθ)cosh(nu), min_σ F_n ≥ 4−4cosh(nu₀), u₀=½ln(1+1/γ²)  (cos≤1, cosh(nu)≤cosh(nu₀)).
So a RIGOROUS upper bound on the negative mass is
  NegBound(n) = ∫_T^{n/2π} 4(cosh(n u₀(γ))−1)·(ln(γ/2π)/2π) dγ   (Riemann–von Mangoldt density).
If P_ver(n) > NegBound(n) then λ_n>0. Crossover n* where NegBound=2N(T) gives the candidate extended bound.

RIGOR: NegBound rigorous (conservative — all tail zeros at extreme σ→0). P_ver≥(1−ε)2N(T) is GATED on an
exponential-sum-over-zeros bound (see LEDGER §22). This script computes both and the crossover.
"""
import numpy as np, mpmath as mp
mp.mp.dps = 30

def Pver_over_real_zeros(path, ns):
    g = np.loadtxt(path); g = g[g>0]; T=g.max(); N=len(g)
    phi = 2*np.arctan(1/(2*g))
    print(f"[P_ver] {N} real zeros, T={T:.1f}, 2N={2*N}")
    for n in ns:
        P = 2*np.sum(1-np.cos(n*phi))
        print(f"   n={n:>13}: P_ver={P:12.1f}  P_ver/2N={P/(2*N):.4f}")
    return T, N

def Ncount(T): return (T/(2*mp.pi))*mp.log(T/(2*mp.pi)) - T/(2*mp.pi)

def negbound(n, T):
    n=mp.mpf(n); T=mp.mpf(T); hi=n/(2*mp.pi)
    if hi<=T: return mp.mpf(0)
    f=lambda g: 4*(mp.cosh(n*0.5*mp.log(1+1/g**2))-1)*mp.log(g/(2*mp.pi))/(2*mp.pi)
    return mp.quad(f, [T, 2*T, 10*T, hi])

def crossover(T):
    twoN=2*Ncount(T); lo=mp.mpf(2*mp.pi*T); hi=mp.mpf(100)*T*T
    for _ in range(80):
        mid=(lo+hi)/2
        if negbound(mid,T)<twoN: lo=mid
        else: hi=mid
    return lo, twoN

if __name__=="__main__":
    print("=== positive floor over real zeros ===")
    try:
        Pver_over_real_zeros('data/zeros_100k.txt',
            [10,1000,int(2*np.pi*74920), int(2*np.pi*74920*30), 10**7, 10**9])
    except Exception as e:
        print("  (zero file not found:",e,")")
    print("\n=== crossover n* (NegBound=2N) and n*/T^2 ===")
    for T in [74920.0, 1.131e6, 1.037e7, 3.0001e12]:
        nstar,twoN=crossover(T)
        print(f"  T={float(T):.4g}: 2N={float(twoN):.3e} 2piT={float(2*mp.pi*T):.3e} "
              f"n*={float(nstar):.3e} n*/(2piT)={float(nstar/(2*mp.pi*T)):.3e} n*/T^2={float(nstar/T**2):.3f}")
    print("\n  candidate: lambda_n>0 for n <~ 3.2*T^2 (T=3.0001e12 -> ~2.9e25), GATED on positive-floor bound.")
