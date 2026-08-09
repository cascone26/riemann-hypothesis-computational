#!/usr/bin/env python3
"""§861: the PRIME-side structure of the RH-critical near-null mode (dual to §858's zero-side). Extract the
margin eigenvector v; decompose its Weil value Q(v) into arch(v), pole(v), and per-prime contributions
<v, M_n v>. Which primes does the RH-critical mode feel most? Compare to §858 (zeros) and §860 (sensitivity)."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
def vonmangoldt(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
NMAX=3000; LAM=vonmangoldt(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW0=LAM[2:]/np.sqrt(NS)
a,b,W=14.0,44.0,1.0   # 4W=4.0 (near-null modes present)
fourW=4*W; N=13; dr=0.004
RGRID=np.arange(-(b+120),(b+120),dr); OMEGA=np.real(digamma(0.25+0.5j*RGRID))-np.log(np.pi)
CEN=np.linspace(a,b,N)
def fej(x):
    y=W*x; out=np.ones_like(y); nz=np.abs(y)>1e-12; out[nz]=(np.sin(y[nz])/y[nz])**2; return out
P=np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G=(P@P.T)*dr
arch=(P*OMEGA[None,:])@P.T*dr/(2*np.pi)
def fc(c):
    y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
mask=LOGN<fourW; ln=LOGN[mask]; nn=NS[mask]; w0=LAMW0[mask]; cosu=np.cos(np.outer(ln,RGRID))
Mn=[]
for jj in range(len(ln)):
    col=cosu[jj]; Mk=np.zeros((N,N))
    for j in range(N): Mk[j,:]=2.0*(((P[j][None,:]*P)@col)*dr/(2*np.pi))
    Mn.append(0.5*(Mk+Mk.T))
prime=np.tensordot(w0,np.array(Mn),axes=1)
w,V=eigh(pole+arch-prime,G)
v=V[:,0]  # margin eigenvector (already G-normalized)
Qv=w[0]
av=float(v@arch@v); pv=float(v@pole@v)
print(f"§861 prime-side structure of RH-critical mode, [14,44] 4W=4, N={N}; margin Q(v)={Qv:+.4e}")
print(f"decomposition: arch(v)={av:+.5f}  pole(v)={pv:+.5f}  prime(v)={av+pv-Qv:+.5f}  (arch+pole-prime=Q ✓ {av+pv-(av+pv-Qv):+.2e})")
print(f"{'n':>4} {'weight':>8} {'per-prime <v,M_n v>':>20} {'% of prime(v)':>13}")
primev=av+pv-Qv
contribs=[(nn[i], w0[i], w0[i]*float(v@Mn[i]@v)) for i in range(len(ln))]
for n_,wt,c in sorted(contribs,key=lambda t:-abs(t[2])):
    print(f"{n_:>4} {wt:>8.4f} {c:>20.5f} {100*c/primev:>12.1f}%")
print("Shows which primes the RH-critical near-null mode feels most (prime-side dual of §858's zero-side).")
