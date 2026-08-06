#!/usr/bin/env python3
"""GRH extension of §772/§780: does the Dirichlet L(s,χ₄) Weil form = Gram matrix of ITS zeros?
χ₄ = odd real character mod 4 (parity a=1, conductor q=4). Weil form on φ_a(x)=e^{−x²/2}cos(ax):
  arch_χ(a,b) = (1/2π)∫ φ̂_a φ̂_b Ω_χ dr,  Ω_χ(r) = log(q/π) + Re ψ(1/4 + a/2 + ir/2) = log(4/π)+Re ψ(3/4+ir/2)
  pole = 0                                   (non-principal χ ⇒ L(s,χ) entire)
  prime_χ(a,b) = 2 Σ_n Λ(n)χ₄(n)/√n · g(a,b,log n)   (χ-twisted von Mangoldt)
Test Q_χ(a,b) =? 2 Σ_{γ>0} φ̂_a(γ)φ̂_b(γ)  over the L(s,χ₄) zeros (Lchi4_zeros_complete.npy).
If it matches, §772's Gram identity (Weil form = Gram of zeros) extends to GRH — and by §780 the form
would likewise detect off-line χ-zeros. Verify-before-ship: a wrong arch constant breaks the match."""
import numpy as np, mpmath as mp
mp.mp.dps=30
PI=np.pi
Zc=np.load('Lchi4_zeros_complete.npy')   # positive ordinates of L(s,χ₄) zeros
def chi4(n):
    r=n%4
    return 1 if r==1 else (-1 if r==3 else 0)
# arch weight for χ₄ (parity 1): Ω(r)=log(4/π)+Re ψ(3/4+ir/2)
R=240.0; NG=48000; rg=np.linspace(-R,R,NG)
print("precomputing Ω_χ(r)...", flush=True)
OM=np.array([float(mp.re(mp.digamma(mp.mpf(3)/4+1j*mp.mpf(r)/2)))+np.log(4/PI) for r in rg])
def phihat(a,r): s=np.sqrt(2*PI); return s/2*(np.exp(-(r-a)**2/2)+np.exp(-(r+a)**2/2))
def arch(aj,ak): return np.trapezoid(phihat(aj,rg)*phihat(ak,rg)*OM,rg)/(2*PI)
from sympy import primerange
NM=300000
# χ-twisted von Mangoldt: n=p^k, Λ(n)=log p, weight χ₄(n)=χ₄(p)^k (completely multiplicative)
_N=[];_L=[];_C=[]
for p in primerange(2,NM):
    pk=p;lp=np.log(p);cp=chi4(p)
    while pk<NM:
        _N.append(pk);_L.append(lp);_C.append(chi4(pk));pk*=p
_N=np.array(_N,float);_L=np.array(_L);_C=np.array(_C,float);_LN=np.log(_N);_IS=1/np.sqrt(_N)
def g(aj,ak,y):
    s=aj+ak;d=aj-ak
    return np.exp(-y**2/4)*(np.sqrt(PI)/2)*(np.exp(-s**2/4)*np.cos(d*y/2)+np.exp(-d**2/4)*np.cos(s*y/2))
def prime(aj,ak): return 2*np.sum(_L*_C*_IS*g(aj,ak,_LN))
def Q_arith(aj,ak): return arch(aj,ak)-prime(aj,ak)   # no pole
def Q_gram(aj,ak):  return 2*np.sum(phihat(aj,Zc)*phihat(ak,Zc))

print(f"{'a':>5} {'Q_arith':>13} {'Q_gram(χ₄ zeros)':>18} {'ratio':>8}")
for a in [6,8,10,12,15,20,25]:
    qa=Q_arith(a,a); qg=Q_gram(a,a)
    print(f"{a:>5} {qa:>13.6f} {qg:>18.6f} {qa/qg if qg else float('nan'):>8.4f}")
# full small matrix: is it PD (GRH-consistent) and does arith=gram?
A=[10,14,18,22,26]
Qa=np.array([[Q_arith(A[i],A[j]) for j in range(len(A))] for i in range(len(A))])
Qg=np.array([[Q_gram(A[i],A[j])  for j in range(len(A))] for i in range(len(A))])
print(f"\n5×5 basis {A}: max|Q_arith-Q_gram|={np.max(np.abs(Qa-Qg)):.3e}")
print(f"  min-eig(Q_arith)={np.linalg.eigvalsh(Qa).min():+.4e}  min-eig(Q_gram)={np.linalg.eigvalsh(Qg).min():+.4e} (GRH⇒≥0)")
