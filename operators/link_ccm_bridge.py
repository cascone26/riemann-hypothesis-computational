#!/usr/bin/env python3
"""§303: FIRST direct comparison of the CCM positive form C (|zeta|^2 measure, ≥0 by construction) and the Weil
form Q (=RH) on the SAME Gaussian band. §262 names the gap as the log-derivative ladder |zeta|^2 -> log -> d/ds -> zeta'/zeta.
This instruments it: build both forms, compare spectra, and check the sign-structure of R = C - Q (where R goes
negative = where CCM positivity FAILS to dominate Weil = the RH locus, made concrete). numpy float (structural probe).
Charisma: REPORT."""
import numpy as np
from numpy.polynomial import polynomial as _p

w_gauss=3.0
def F(c,r): return np.exp(-(r-c)**2/(2*w_gauss**2))+np.exp(-(r+c)**2/(2*w_gauss**2))
def Omega(r):
    from scipy.special import digamma
    return np.real(digamma(0.25+1j*r/2))-np.log(np.pi)
# von Mangoldt for Weil
def mangoldt(N):
    L=np.zeros(N+1)
    for p in range(2,N+1):
        if all(p%d for d in range(2,int(p**0.5)+1)):
            q=p
            while q<=N: L[q]=np.log(p); q*=p
    return L
lo,hi,M=14.0,38.0,9
centers=np.linspace(lo,hi,M)
# r-grid for all integrals
rg=np.linspace(0.0,60.0,6001); dr=rg[1]-rg[0]
Fmat=np.array([F(c,rg) for c in centers])            # M x len(rg)
# --- Weil form Q = pole + arch - prime ---
Om=Omega(rg)
Np=60; L=mangoldt(Np)
def trap(y): return np.trapezoid(y,rg)
Q=np.zeros((M,M)); C=np.zeros((M,M))
# CCM positive weight w_ccm(r)=|sum_{n<=M0} n^{-1/2} e^{-i r ln n}|^2  (truncated |zeta(1/2+ir)|^2, >=0)
M0=300
ns=np.arange(1,M0+1); coef=ns**-0.5; logn=np.log(ns)
Zt=(coef[:,None]*np.exp(-1j*np.outer(logn,rg))).sum(axis=0)   # truncated zeta on line
w_ccm=np.abs(Zt)**2
# normalize CCM weight to comparable scale (unit mean over band-support) so C,Q are comparable in magnitude
supp=(rg>lo-12)&(rg<hi+12)
w_ccm_n=w_ccm/np.mean(w_ccm[supp])
for j in range(M):
    for k in range(M):
        fk=Fmat[j]*Fmat[k]
        arch=(1/(2*np.pi))*2*trap(fk*Om)
        pole=2*( (np.exp(-(0.5j-centers[j])**2/(2*w_gauss**2))+np.exp(-(0.5j+centers[j])**2/(2*w_gauss**2))) *
                 (np.exp(-(0.5j-centers[k])**2/(2*w_gauss**2))+np.exp(-(0.5j+centers[k])**2/(2*w_gauss**2))) ).real
        prime=2*sum(L[n]/np.sqrt(n)*(1/(2*np.pi))*2*trap(fk*np.cos(np.log(n)*rg)) for n in range(2,Np+1) if L[n]>0)
        Q[j,k]=pole+arch-prime
        C[j,k]=(1/(2*np.pi))*2*trap(fk*w_ccm_n)     # CCM positive form (|zeta|^2 weight)
Q=(Q+Q.T)/2; C=(C+C.T)/2
eQ=np.linalg.eigvalsh(Q); eC=np.linalg.eigvalsh(C)
print("Weil Q eigs:",np.array2string(eQ,precision=4,suppress_small=True),flush=True)
print("CCM  C eigs:",np.array2string(eC,precision=4,suppress_small=True),"(all>=0 => positive form, as expected)",flush=True)
print(f"  Q min={eQ.min():.3e} (>=0 => RH-in-band)  C min={eC.min():.3e} (>=0 by construction)",flush=True)
# R = C - Q : if R <=0 (neg semidef) everywhere => C dominates Q => Q<=C (upper bd, not RH).
# if R >=0 => Q<=C... we want sign structure. The RH-relevant question: is Q >= 0 recoverable from C>=0?
for scaleC in [1.0]:
    R=scaleC*C-Q; eR=np.linalg.eigvalsh(R)
    print(f"  R=C-Q eigs:",np.array2string(eR,precision=4,suppress_small=True),flush=True)
    print(f"    R indefinite? min={eR.min():.3e} max={eR.max():.3e} -> {'INDEFINITE (linear bridge fails, as §262 predicts)' if eR.min()<-1e-6 and eR.max()>1e-6 else 'DEFINITE'}",flush=True)
# The REAL bridge (§262): Weil weight = (log n /2) x (log-CCM spectral weight). Verify the identity on prime powers:
print("bridge identity check  Lambda(n)/sqrt(n)  vs  (log n/2)*[log|zeta|^2 coeff (2/k)p^-k/2]:",flush=True)
for (n,) in [(2,),(3,),(4,),(5,),(8,),(9,)]:
    # n=p^k
    p=2
    for pp in [2,3,5,7]:
        kk=1; q=pp
        while q<n: q*=pp; kk+=1
        if q==n: p=pp; k=kk; break
    else:
        continue
    weil=L[n]/np.sqrt(n)
    logccm=(2.0/k)*p**(-k/2.0)   # coeff of log|zeta|^2 at freq k log p
    bridge=(np.log(n)/2.0)*logccm
    print(f"  n={n}: Weil={weil:.6f}  (logn/2)*logCCM={bridge:.6f}  ratio={weil/bridge:.6f}",flush=True)
print("DONE",flush=True)
