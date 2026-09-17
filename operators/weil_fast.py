import numpy as np
import mpmath as mp
def digamma(z): return complex(mp.digamma(z))
def trigamma(z): return complex(mp.polygamma(1,z))
from scipy.special import digamma as _rdig
def prime_powers_up_to(c):
    primes=[];out=[]
    x=2
    while x<=c:
        if all(x%p for p in primes): primes.append(x)
        x+=1
    for p in primes:
        q=p
        while q<=c: out.append((q,p)); q*=p
    return out
def build_fast(c,N):
    L=np.log(c); PI=np.pi; sp2=16*PI*PI; l2=L*L
    pref02=32*L*np.sinh(L/4)**2
    # arch closed forms (float): S,CC,XC via digamma/trigamma + geometric tail
    ns=np.arange(0,N+1); 
    def geom(n):
        w=2*PI*n/L; w2=w*w; gS=gCC=gX1=gX2=0.0
        for k in range(0,80):
            ck=2*k+0.5; e=np.exp(-ck*L); den=ck*ck+w2
            gS+=e/den
            if n!=0: gCC+=e*w2/(ck*den)
            gX1+=e*ck/den; gX2+=e*(ck*ck-w2)/(den*den)
        return gS,gCC,gX1,gX2
    S=np.zeros(N+1);CC=np.zeros(N+1);XC=np.zeros(N+1)
    psiq=float(mp.digamma(mp.mpf(1)/4))
    for n in ns:
        w=2*PI*n/L; z=0.25+1j*PI*n/L
        psi=digamma(z); psi1=trigamma(z)
        gS,gCC,gX1,gX2=geom(n)
        if n==0: S[n]=0;CC[n]=0
        else:
            S[n]=0.5*psi.imag - w*gS
            CC[n]=-0.5*(psi.real-psiq)+gCC
        XC[n]=0.25*psi1.real - L*gX1 - gX2
    kappa=np.log(4*PI*(np.exp(L)-1)/(np.exp(L)+1))+np.euler_gamma
    U=np.exp(L/2); J=-2*np.log(U+1)+np.log(U*U+1)+2*np.arctan(U)+np.log(2)-PI/2
    # prime side: precompute A_k, Bdiag_k
    pdata=prime_powers_up_to(c)
    y=np.array([np.log(q) for q,p in pdata]); w=np.array([np.log(p)/np.sqrt(q) for q,p in pdata])
    ks=np.arange(-N,N+1)
    A=np.array([np.sum(w*np.sin(2*PI*k*y/L)) for k in ks])          # A_k
    Bd=np.array([np.sum(w*2*(1-y/L)*np.cos(2*PI*k*y/L)) for k in ks])# diag prime
    DIM=2*N+1; T=np.zeros((DIM,DIM))
    idx=lambda n: n+N
    for i in range(DIM):
        n=i-N
        for j in range(i,DIM):
            m=j-N
            num=l2-sp2*m*n; den=(l2+sp2*m*m)*(l2+sp2*n*n); W02=pref02*num/den
            if n==m:
                WR=kappa+2*CC[abs(n)]+J-(2/L)*XC[abs(n)]; Wp=Bd[idx(n)]
            else:
                Ss=lambda k: S[k] if k>=0 else -S[-k]
                WR=(Ss(m)-Ss(n))/(PI*(n-m))
                Wp=(A[idx(m)]-A[idx(n)])/(PI*(n-m))
            v=W02-WR-Wp; T[i,j]=v; T[j,i]=v
    return T
if __name__=="__main__":
    import sys
    sys.path.insert(0, '../external/groskin_2607.02828/code')
    from flint import ctx
    from dh_certified_tau import build_arb_tau
    # validate float-fast vs Arb (float) at small c,N
    c,N=25,20; ctx.prec=400
    A,DIM=build_arb_tau(c,N,400); Ma=np.array([[float(A[i,j].mid()) for j in range(DIM)] for i in range(DIM)])
    Mf=build_fast(c,N)
    print(f"c={c},N={N}: max|fast-arb| = {np.max(np.abs(Mf-Ma)):.2e}  (should be ~1e-10)")
