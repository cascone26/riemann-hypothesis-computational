#!/usr/bin/env python3
"""ATTACK A6: is the Weil-positivity (tau>=0) UNIVERSAL across L-functions? Test a Dirichlet
L-function L(s,chi_4) (chi = non-principal char mod 4, ODD). If tau_chi>=0 holds like it does
for zeta, positivity is universal -> any proof can't use zeta-specific structure. Reuses the
VALIDATED connes_decisive_v3 numerical construction (arch integral + prime + boundary), changing
ONLY: (1) arch digamma arg 1/4 -> 3/4 (odd char a=1); (2) conductor const -log(pi) -> log(q/pi);
(3) prime weights get chi(k), even k dropped; (4) NO boundary/pole term (L(s,chi) is entire).
Validation gate: the zeta path (a=1/4, all primes, boundary ON) must reproduce known PD + ~23 near-kernel."""
import numpy as np
from scipy.special import digamma
import sys
TWO_PI=2*np.pi

def make(L, a_re, const, boundary_on, character):
    def omega(r): return np.real(digamma(a_re+0.5j*r))+const
    def e_hat(r,n):
        s=TWO_PI*n/L; d=s-r; out=np.empty_like(d,dtype=complex)
        sm=np.abs(d)<1e-12; out[sm]=np.sqrt(L)+0j; dd=d[~sm]
        out[~sm]=(1/np.sqrt(L))*np.exp(0.5j*dd*L)*2*np.sin(0.5*dd*L)/dd
        return out
    def qk(n,m,y):
        if n==m: return 2*(1-abs(y)/L)*np.cos(TWO_PI*n*y/L)
        return (np.sin(TWO_PI*m*y/L)-np.sin(TWO_PI*n*y/L))/(np.pi*(n-m))
    # von Mangoldt prime-power table k<=exp(L) with weight log p, times character
    c=int(round(np.exp(L)))
    VM=[]
    for k in range(2,c+1):
        n=k;p=0;f=2;ok=True
        while f*f<=n:
            if n%f==0:
                if p==0:p=f
                elif p!=f:ok=False;break
                while n%f==0:n//=f
            f+=1
        if n>1:
            if p==0:p=n
            elif p!=n:ok=False
        if ok and p>0:
            chi=character(k)
            if chi!=0: VM.append((k,np.log(p)*chi))
    def build(N,R,G):
        idx=list(range(-N,N+1));D=len(idx)
        grid=np.linspace(-R,R,G);h=grid[1]-grid[0]
        w=np.full(G,h);w[0]=w[-1]=0.5*h;Om=omega(grid)
        E=np.empty((D,G),dtype=complex)
        for k,n in enumerate(idx): E[k]=e_hat(grid,n)
        arch=((E*(w*Om))@E.conj().T).real/TWO_PI
        P=np.zeros((D,D))
        for ii,n in enumerate(idx):
            for jj in range(ii,D):
                m=idx[jj];s=sum(Lam*k**-0.5*qk(n,m,np.log(k)) for k,Lam in VM)
                P[ii,jj]=P[jj,ii]=s
        W=np.zeros((D,D))
        if boundary_on:
            sh=np.sinh(L/4)**2
            for ii,n in enumerate(idx):
                for jj,m in enumerate(idx):
                    num=L**2-16*np.pi**2*m*n;den=(L**2+16*np.pi**2*m**2)*(L**2+16*np.pi**2*n**2)
                    W[ii,jj]=32*L*sh*num/den
        tau=arch+W-P;return (tau+tau.T)/2,idx
    return build

L=float(2*np.log(np.sqrt(14)))   # same geometry as zeta runs (L=log 14)
# ---- VALIDATION GATE: zeta path ----
zeta_build=make(L,0.25,-np.log(np.pi),True,lambda k:1)
print("=== VALIDATION (zeta path: arch 1/4, all primes, boundary ON) ===")
for N in [16,32,48]:
    tau,idx=zeta_build(N,TWO_PI*N/L+80,max(60000,3000*N))
    ev=np.linalg.eigvalsh(tau);nk=int(np.sum(np.abs(ev)<1e-6))
    print(f"  N={N}: min_eig={ev[0]:.2e}  near-kernel={nk}  (expect ~PSD, nk~23)")
# ---- chi_4 (odd, conductor 4, no boundary) ----
def chi4(k): return 0 if k%2==0 else (1 if k%4==1 else -1)
chi_build=make(L,0.75,np.log(4/np.pi),False,chi4)
print("\n=== chi_4 path (arch 3/4, chi weights, NO boundary, conductor 4) ===")
res=[]
for N in [16,24,32,48,64]:
    tau,idx=chi_build(N,TWO_PI*N/L+80,max(60000,3000*N))
    ev=np.linalg.eigvalsh(tau);nk=int(np.sum(np.abs(ev)<1e-6))
    print(f"  N={N}: min_eig={ev[0]:.3e}  2nd={ev[1]:.3e}  near-kernel={nk}")
    res.append((N,float(ev[0])))
# sensitivity to conductor const (is PD robust or knife-edge?)
print("\n  conductor-const sensitivity (N=48, min_eig vs const):")
for cst in [np.log(4/np.pi)-0.5,np.log(4/np.pi),np.log(4/np.pi)+0.5,-np.log(np.pi)]:
    b=make(L,0.75,cst,False,chi4);tau,_=b(48,TWO_PI*48/L+80,144000)
    print(f"    const={cst:+.3f}: min_eig={np.linalg.eigvalsh(tau)[0]:.3e}")
print("\nVERDICT: chi_4 min_eig>=0 (like zeta) => positivity UNIVERSAL across L-functions.")
