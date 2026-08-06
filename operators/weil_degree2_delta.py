#!/usr/bin/env python3
"""DEGREE-2 (GL(2)) test of the Gram framework: does the Weil form of the Ramanujan Δ L-function
(weight 12, level 1) = Gram matrix of ITS zeros? First test beyond degree-1 (ζ/Dirichlet).

Analytic normalization (functional equation s↔1−s): a_p = τ(p)/p^{11/2} (|a_p|≤2, Deligne).
Arch factor Γ_ℂ(s+11/2) ⇒ arch density Ω_Δ(r) = Re ψ(6 + ir) − log(2π)   [s=½+ir, +11/2 → 6].
Prime term (degree 2): Σ_{p^k} b_{p^k}/p^{k/2}·log p · g(k log p), b = Satake power sums
  b_{p^0}=2, b_{p^1}=a_p, b_{p^{k}}=a_p·b_{p^{k-1}} − b_{p^{k-2}}   (α+β=a_p, αβ=1).
Test Q_arith(a,b)=arch−prime  =?  2 Σ_{γ>0} φ̂_a(γ)φ̂_b(γ)  over Δ zeros (Delta_zeros_complete.npy).
Verify-before-ship: a wrong arch shift / Satake convention breaks the match."""
import numpy as np, mpmath as mp
mp.mp.dps=30
PI=np.pi
Zc=np.load('Delta_zeros_complete.npy')   # 228 zeros of L(s,Δ), analytic normalization

# --- τ(n) via Δ/q = ∏(1−qⁿ)²⁴  (truncate at NMAX) ---
NMAX=4000
P=np.zeros(NMAX+1); P[0]=1.0     # ∏(1−qⁿ)
for n in range(1,NMAX+1):
    # multiply P by (1 − q^n): P[k] -= P[k-n]
    P[n:] -= P[:NMAX+1-n]
# raise to 24th power (mod q^{NMAX+1}) by repeated multiply
def pmul(a,b):
    c=np.zeros(NMAX+1)
    for k in range(NMAX+1):
        if a[k]!=0:
            c[k:]+=a[k]*b[:NMAX+1-k]
    return c
Q=np.array([1.0]+[0]*NMAX)
base=P.copy(); e=24
while e:
    if e&1: Q=pmul(Q,base)
    e>>=1
    if e: base=pmul(base,base)
# Δ = q·Q  ⇒ τ(n) = Q[n-1]
tau=np.zeros(NMAX+1)
tau[1:]=Q[:NMAX]
print(f"τ(1..7) = {tau[1:8]}  (expect 1,-24,252,-1472,4830,-6048,-16744)", flush=True)

def ap(p): return tau[p]/p**5.5   # normalized Hecke eigenvalue, |ap|<=2
from sympy import primerange
# Satake von Mangoldt data: for each p^k<=NMAX, weight b_{p^k}=α^k+β^k, freq log(p^k)
_L=[];_W=[]
for p in primerange(2,NMAX):
    a=ap(p); b0,b1=2.0,a; k=1; pk=p
    while pk<=NMAX:
        bk = b1 if k==1 else (a*b1-b0)  # b_{p^k}
        if k>=2: b0,b1=b1,bk
        _L.append(np.log(pk)); _W.append((bk if k>=1 else 2.0)*np.log(p)/np.sqrt(pk))
        # advance
        if k==1: pass
        k+=1; pk*=p
        if k==2: b0,b1=a,(a*a-2)  # reset recursion cleanly below
# recompute cleanly (the inline recursion above is fiddly) — do it straight:
_L=[];_W=[]
for p in primerange(2,NMAX):
    a=ap(p); lp=np.log(p)
    b=[2.0,a]  # b[0]=b_{p^0}, b[1]=b_{p^1}
    k=1; pk=p
    while pk<=NMAX:
        if k>=2: b.append(a*b[k-1]-b[k-2])
        _L.append(np.log(pk)); _W.append(b[k]*lp/np.sqrt(pk))
        k+=1; pk*=p
_L=np.array(_L);_W=np.array(_W)

R=240.0;NG=48000;rg=np.linspace(-R,R,NG)
print("precomputing Ω_Δ(r)=Re ψ(6+ir)−log(2π)...", flush=True)
OM=np.array([float(mp.re(mp.digamma(6+1j*mp.mpf(r))))-np.log(2*PI) for r in rg])
def ph(a,r):s=np.sqrt(2*PI);return s/2*(np.exp(-(r-a)**2/2)+np.exp(-(r+a)**2/2))
def arch(i,j):return np.trapezoid(ph(i,rg)*ph(j,rg)*OM,rg)/(2*PI)
def g(i,j,y):
    s=i+j;d=i-j
    return np.exp(-y**2/4)*(np.sqrt(PI)/2)*(np.exp(-s**2/4)*np.cos(d*y/2)+np.exp(-d**2/4)*np.cos(s*y/2))
def prime(i,j): return np.sum(_W*g(i,j,_L))   # note: NOT ×2 (degree-2 sum already over p^k with b)
def Qa(i,j): return 2*(arch(i,j)-prime(i,j))   # ×2 = Γ_ℂ degree-2 normalization (Γ_ℂ=Γ_ℝ·Γ_ℝ)
def Qg(i,j): return 2*np.sum(ph(i,Zc)*ph(j,Zc))
print('\nL(s,Δ) degree-2 GRH: Q_arith vs Q_gram(Δ zeros)')
for a in [10,12,14,16,18,20]:
    print(f'  a={a}: arith={Qa(a,a):+.5f} gram={Qg(a,a):+.5f} ratio={Qa(a,a)/Qg(a,a):.4f}')
