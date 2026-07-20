#!/usr/bin/env python3
"""ATTACK A13: is the SUSY superpotential of tau^{1/2} CLEAN in the archimedean (A+B) eigenbasis?
A12 found tau^{1/2} quasi-local w/ jagged superpotential in the frequency basis. Test whether a
natural basis (eigenbasis of arch+boundary A+B, where tau = Diag(gaps) - P_coupling) smooths W."""
import numpy as np
from scipy.special import digamma
import mpmath as _mp
def poly1(z): return complex(_mp.polygamma(1,complex(z)))

def build(N,c,prime_on=True):
    L=np.log(c);PI=np.pi;psi_q=digamma(0.25)
    S=np.zeros(N+1);CC=np.zeros(N+1);XC=np.zeros(N+1)
    for n in range(N+1):
        w=2*PI*n/L;z=0.25+1j*PI*n/L;psi=digamma(z);psi1=poly1(z)
        ks=np.arange(0,80);ck=2*ks+0.5;e=np.exp(-ck*L);den=ck*ck+w*w
        gS=np.sum(e/den);gCC=np.sum(e*w*w/(ck*den)) if n else 0.0
        gX1=np.sum(e*ck/den);gX2=np.sum(e*(ck*ck-w*w)/(den*den))
        if n: S[n]=psi.imag/2-w*gS;CC[n]=-(psi.real-psi_q)/2+gCC
        XC[n]=psi1.real/4-L*gX1-gX2
    sp2=16*PI*PI;l2=L*L;pref=32*L*np.sinh(L/4)**2
    eL=np.exp(L);kap=np.log(4*PI*(eL-1)/(eL+1))+np.euler_gamma
    U=np.exp(L/2);J=-2*np.log(U+1)+np.log(U*U+1)+2*np.arctan(U)+np.log(2)-PI/2
    if prime_on:
        pr=[];x=2
        while x<=c:
            if all(x%p for p in range(2,int(x**.5)+1)):
                q=x
                while q<=c: pr.append((q,x));q*=x
            x+=1
        wts=np.array([np.log(p)*q**-.5 for q,p in pr]);pos=np.array([np.log(q) for q,p in pr])
    else: wts=np.array([]);pos=np.array([])
    def Sg(nn): return S[nn] if nn>=0 else -S[-nn]
    D=2*N+1;T=np.zeros((D,D));ns=np.arange(-N,N+1)
    for i in range(D):
        n=ns[i]
        for j in range(i,D):
            m=ns[j];num=l2-sp2*m*n;dd=(l2+sp2*m*m)*(l2+sp2*n*n);W02=pref*num/dd
            WR=kap+2*CC[abs(n)]+J-(2/L)*XC[abs(n)] if n==m else (Sg(m)-Sg(n))/(PI*(n-m))
            Wp=0.0
            if len(wts):
                qv=2*(1-pos/L)*np.cos(2*PI*n*pos/L) if n==m else (np.sin(2*PI*m*pos/L)-np.sin(2*PI*n*pos/L))/(PI*(n-m))
                Wp=np.sum(wts*qv)
            v=W02-WR-Wp;T[i,j]=v;T[j,i]=v
    return T

def sq(M): w,V=np.linalg.eigh(M);return V@np.diag(np.sqrt(np.maximum(w,0)))@V.T
def band(M):
    D=len(M);return np.array([np.mean([abs(M[n,n+d]) for n in range(D-d)]) for d in range(min(D,20))])

N,c=64,14
tau=build(N,c,True); AB=build(N,c,False)
wAB,Q=np.linalg.eigh(AB)                 # arch+boundary eigenbasis
tp=Q.T@tau@Q                             # tau in AB-eigenbasis
St_freq=sq(tau); St_AB=sq(tp)
bf=band(np.abs(St_freq));bf/=bf[0]; ba=band(np.abs(St_AB));ba/=ba[0]
# superpotential smoothness = variation of consecutive diagonal (2nd difference RMS, bulk only)
def smooth(M):
    d=np.abs(np.diag(M));d=d[d>1e-4]  # drop near-kernel zeros
    return float(np.sqrt(np.mean(np.diff(d,2)**2))/np.mean(d))
print("=== A13: superpotential in frequency basis vs A+B-eigenbasis (N=64,c=14) ===")
print(f" AB min eig (arch+boundary indefinite? proj said 4 neg): {wAB.min():.3f}  #neg={int(np.sum(wAB<0))}")
print(f" tau^1/2 offdiag decay  FREQ basis: {[f'{x:.3f}' for x in bf[:7]]}")
print(f" tau^1/2 offdiag decay  A+B  basis: {[f'{x:.3f}' for x in ba[:7]]}")
print(f" superpotential jaggedness (2nd-diff RMS/mean, lower=cleaner):")
print(f"   FREQ basis: {smooth(St_freq):.3f}    A+B basis: {smooth(St_AB):.3f}")
print(f" band-energy within |d|<=1  FREQ: {np.linalg.norm(np.triu(St_freq,-1)-np.triu(St_freq,2))/np.linalg.norm(St_freq):.3f}")
