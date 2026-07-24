#!/usr/bin/env python3
"""
cvs_tau_operator.py — LEDGER §165 (frontier-grade tool rebuild): the EXACT Connes–van Suijlekom truncated Weil matrix,
targeting Groskin's published λ_min^even(c=13,N=100)=2.865e-59 as the correctness check.
CCM 2511.22755 Prop 4.3 / Lemma 5.1 + Groskin Eq 5:
  b_n = −(1/π) ∫₀^L sin(2πnx/L) ρ(x) dx,   ρ(x)=e^{x/2}/(e^x−e^{−x}),   L=log c   (b_{-n}=−b_n, b_0=0)
  a_n = b'(n) = −(2/L) ∫₀^L x cos(2πnx/L) ρ(x) dx   (diagonal = divided-difference limit; a_{-n}=a_n)
  τ_{ij} = a_i (i=j) else (b_i−b_j)/(i−j),  indices n=−N..N (dim 2N+1).
Test: smallest eigenvalue vs N (does it PLUMMET → toward 2.865e-59 at N=100?). ρ prime-free ⟹ if it plummets, ρ is the
complete CCM density (arithmetic absorbed in the construction); if it gives only an arch model (no plummet), add the
prime part b_n^prime = −(1/π) Σ_{p^m≤c}(log p)p^{−m/2} sin(2πn·m·log p/L) (CCM Eq 3.16). High precision (mpmath).
"""
import sys, mpmath as mp

def build_and_smallest(c, N, dps, add_primes=False):
    mp.mp.dps = dps
    L = mp.log(c)
    def rho(x):
        if x < mp.mpf(10)**(-(dps//2)): return mp.mpf(0)   # integrand sin*rho→finite; guard the 1/x pole (sin kills it)
        return mp.e**(x/2)/(mp.e**x - mp.e**(-x))
    def bn(n):
        f = lambda x: mp.sin(2*mp.pi*n*x/L)*rho(x)
        val = -(1/mp.pi)*mp.quad(f, [0, L/8, L/4, L/2, L])
        if add_primes:
            # prime part (CCM Eq 3.16): sum over prime powers p^m <= c
            s = mp.mpf(0); P = int(c)
            sieve=[True]*(P+1)
            for p in range(2,P+1):
                if sieve[p]:
                    for q in range(p*p,P+1,p): sieve[q]=False
                    m=1; pk=p
                    while pk<=P:
                        s += mp.log(p)*mp.mpf(pk)**(mp.mpf(-1)/2)*mp.sin(2*mp.pi*n*(m*mp.log(p))/L)
                        m+=1; pk*=p
            val += -(1/mp.pi)*s
        return val
    def an(n):
        f = lambda x: x*mp.cos(2*mp.pi*n*x/L)*rho(x)
        return -(2/L)*mp.quad(f, [0, L/8, L/4, L/2, L])
    # precompute b_n, a_n for n=0..N
    b = [mp.mpf(0)] + [bn(n) for n in range(1,N+1)]
    a = [an(n) for n in range(0,N+1)]
    idx = list(range(-N, N+1))
    def B(n): return -b[-n] if n<0 else b[n]
    def A(n): return a[abs(n)]
    T = mp.zeros(2*N+1)
    for ii,i in enumerate(idx):
        for jj,j in enumerate(idx):
            T[ii,jj] = A(i) if i==j else (B(i)-B(j))/(i-j)
    T = (T+T.T)/2
    E = mp.eigsy(T, eigvals_only=True)
    return [E[k] for k in range(min(6,len(E)))]

c=13; dps=int(sys.argv[1]) if len(sys.argv)>1 else 60
addp = (len(sys.argv)>2 and sys.argv[2]=='primes')
print(f"§165 CvS τ-operator, c={c} (L=log13={float(mp.log(13)):.4f}), dps={dps}, primes={addp}")
print(f"target: Groskin λ_min^even(c=13,N=100)=2.865e-59\n")
print(f"{'N':>5} {'smallest eigenvalues (lowest 4)':>60}")
for N in [12, 20, 30]:
    ev = build_and_smallest(c, N, dps, add_primes=addp)
    print(f"{N:>5}  " + "  ".join(mp.nstr(e,4) for e in ev[:4]), flush=True)
print("\nREAD: if the smallest |eigenvalue| PLUMMETS as N grows (12→20→30, dropping many orders) → the τ-operator is right")
print("and ρ (Prop 4.3) is the complete density; scale to N=100 (Olympus) & check vs 2.865e-59. If it stays O(1)/doesn't")
print("plummet → ρ is arch-only; rerun with 'primes' to add the CCM Eq 3.16 prime part.")
