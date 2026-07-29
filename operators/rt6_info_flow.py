#!/usr/bin/env python3
"""§440 (prime↔zero information flow, novel): how many zeros K are needed to resolve prime p to within tol? Recover
Λ(p) using only the first K zeros (§438 method), find K_min(p). Prediction (Fourier/uncertainty): need zeros up to
γ~p, i.e. K_min ≈ N(p) (Riemann-von Mangoldt count up to height p). Quantifies how the arithmetic is distributed
across the zero sequence."""
import numpy as np
from scipy.special import digamma
from scipy import integrate
Zall=np.load('zeros150.npy')
def Om(r): return np.real(digamma(0.25+0.5j*r))-np.log(np.pi)
def recover(n0, Z):
    L=np.log(n0); gap=min(L-np.log(n0-1),np.log(n0+1)-L); sig=0.40*gap
    h=lambda r: np.exp(-sig**2*r**2/2)*np.cos(r*L)
    lhs=2*np.sum(h(Z))
    arch,_=integrate.quad(lambda r: h(r)*Om(r),-600,600,limit=900); arch/=(2*np.pi)
    pole=2*np.exp(sig**2/8)*np.cosh(L/2); P=(arch+pole-lhs)/2
    return P/(1/(2*np.sqrt(2*np.pi)*sig))*np.sqrt(n0)
# calibrate const on full data n=2
cal=np.log(2)/recover(2,Zall)
def Nrvm(T): 
    f=lambda t:(t/(2*np.pi))*np.log(t/(2*np.pi))-t/(2*np.pi); return f(T)
print("§440 K_min(p): #zeros to recover Λ(p) within 10%. predict K_min≈N(p) (RvM count to height p).")
print(f"{'p':>4} {'K_min':>6} {'γ_Kmin':>8} {'p (target γ)':>12} {'N(p) RvM':>9}")
for p in [3,5,7,11,13,17,23,31,43]:
    Kmin=None
    for K in range(5,len(Zall)+1,2):
        rec=recover(p,Zall[:K])*cal
        if abs(rec-np.log(p))/np.log(p) < 0.10:
            # require stability: also next K within tol
            rec2=recover(p,Zall[:min(K+4,len(Zall))])*cal
            if abs(rec2-np.log(p))/np.log(p)<0.12: Kmin=K; break
    if Kmin:
        print(f"{p:>4} {Kmin:>6} {Zall[Kmin-1]:8.2f} {p:12d} {Nrvm(p):9.1f}",flush=True)
    else:
        print(f"{p:>4} {'>150':>6} {'-':>8} {p:12d} {Nrvm(p):9.1f}  (not resolved w/ 150 zeros)")
print("""
reading: if K_min tracks N(p) (RvM count to height p) i.e. γ_Kmin ≈ p, the arithmetic of prime p is encoded in
the zeros UP TO height ~p — a clean information-flow law (uncertainty principle: resolving frequency log p needs
the zeros as a 'signal' out to γ~p). Confirms the local, quantitative prime↔zero correspondence.
""")
