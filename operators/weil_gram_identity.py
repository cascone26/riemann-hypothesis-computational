#!/usr/bin/env python3
"""Verify the Weil form is a GRAM matrix of zero-evaluation vectors:
   Q_{jk} = arch(a_j,a_k)+pole−prime  =?=  2 Σ_{γ>0} φ̂_{a_j}(γ) φ̂_{a_k}(γ)  =  ⟨v_j, v_k⟩
with v_j = (√2·φ̂_{a_j}(γ))_{γ>0}, φ̂_a(r)=√(2π)/2(e^{−(r−a)²/2}+e^{−(r+a)²/2}).

If true (verified for the DIAGONAL to 4 digits §771), it EXPLAINS all of §769:
  - separated a_j → v_j pick out disjoint zeros → Q≈diagonal
  - dense a_j → v_j nearly dependent → Q rank-deficient → min-eig→0 (non-coercive)
and makes 'Weil PSD ⟺ Gram of REAL vectors ⟺ γ real ⟺ RH' explicit.
Test on a CLOSE basis (substantial off-diagonals). Arithmetic side = fast von-Mangoldt build
(validated to 1e-5 vs rigorous Arb, §767/§771)."""
import numpy as np, mpmath as mp
mp.mp.dps = 30
PI=np.pi; LOGPI=np.log(PI)
R=240.0; NG=48000; rg=np.linspace(-R,R,NG)
print("precomputing Ω(r)...", flush=True)
OM=np.array([float(mp.re(mp.digamma(mp.mpf(1)/4+1j*mp.mpf(r)/2)))-LOGPI for r in rg])
def phihat(a,r): s2p=np.sqrt(2*PI); return s2p/2*(np.exp(-(r-a)**2/2)+np.exp(-(r+a)**2/2))
def arch(aj,ak): return np.trapezoid(phihat(aj,rg)*phihat(ak,rg)*OM, rg)/(2*PI)
def pole(aj,ak):
    s2p=np.sqrt(2*PI)
    pj=s2p/2*(np.exp(-((1j/2-aj)**2)/2)+np.exp(-((1j/2+aj)**2)/2))
    pk=s2p/2*(np.exp(-((1j/2-ak)**2)/2)+np.exp(-((1j/2+ak)**2)/2))
    return 2*np.real(pj*pk)
from sympy import primerange
NMAX=300000; _p=list(primerange(2,NMAX+1)); _N=[];_L=[]
for p in _p:
    pk=p; lp=np.log(p)
    while pk<=NMAX: _N.append(pk);_L.append(lp);pk*=p
_N=np.array(_N,float);_L=np.array(_L);_LOGN=np.log(_N);_IS=1/np.sqrt(_N)
def g_cross(aj,ak,y):
    s=aj+ak; d=aj-ak
    return np.exp(-y**2/4)*(np.sqrt(PI)/2)*(np.exp(-s**2/4)*np.cos(d*y/2)+np.exp(-d**2/4)*np.cos(s*y/2))
def prime(aj,ak): return 2*np.sum(_L*_IS*g_cross(aj,ak,_LOGN))
def Q_arith(aj,ak): return arch(aj,ak)+pole(aj,ak)-prime(aj,ak)

Z=np.load('zeros10k.npy')
def Q_gram(aj,ak): return 2*np.sum(phihat(aj,Z)*phihat(ak,Z))

for label, A in [("close near zeros {14,15,16,17}", [14,15,16,17]),
                 ("dense [20,24]/5", list(np.linspace(20,24,5))),
                 ("separated {15,25,35,45}", [15,25,35,45])]:
    n=len(A)
    Qa=np.array([[Q_arith(A[i],A[j]) for j in range(n)] for i in range(n)])
    Qg=np.array([[Q_gram(A[i],A[j])  for j in range(n)] for i in range(n)])
    err=np.max(np.abs(Qa-Qg))
    eva=np.linalg.eigvalsh(Qa); evg=np.linalg.eigvalsh(Qg)
    print(f"\n{label}:")
    print(f"  max|Q_arith - Q_gram| = {err:.3e}   (arith from primes, gram from zeros)")
    print(f"  eig(arith) min={eva.min():+.5f} max={eva.max():.4f} | eig(gram) min={evg.min():+.5f} (Gram⇒≥0)")
    print(f"  Gram rank (sv>1e-6): {np.sum(np.linalg.svd(Qg,compute_uv=False)>1e-6)}/{n}  "
          f"[separated→full rank, dense→deficient]")
