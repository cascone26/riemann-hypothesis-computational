import numpy as np
from scipy.special import loggamma
def w0(s): return np.exp(2.0*np.real(loggamma(0.25+0.5j*s)))
def wp(s,p): return w0(s)/(1.0-2.0*p**-0.5*np.cos(s*np.log(p))+1.0/p)
def jacobi(weight,L,M,K):
    s=np.linspace(-L,L,M); ds=s[1]-s[0]; w=weight(s)*ds
    Mt=w.sum(); qp=np.zeros_like(s); q=np.sqrt(w/Mt); a=np.zeros(K); b=np.zeros(K); bb=0.0
    for k in range(K):
        v=s*q; aa=np.dot(q,v); a[k]=aa; v=v-aa*q-bb*qp
        v=v-q*np.dot(q,v)-qp*np.dot(qp,v); bb=np.linalg.norm(v)
        if k+1<K: b[k+1]=bb
        qp=q; q=v/bb
    return a,b
def Wmat(a,b,K,lam=np.sqrt(2.0)):
    J=np.diag(a[:K])+np.diag(b[1:K],1)+np.diag(b[1:K],-1)   # mult-by-s (real symmetric Jacobi)
    N=np.diag(np.arange(K))
    return -J@J + 2*np.pi*lam*lam*(4*N+np.eye(K)) - 0.25*np.eye(K)
L,M=50,240000
zeros=np.array([14.1347,21.0220,25.0109,30.4249,32.9351])
print("Absolute low-lying gamma=2*sqrt(|neg eig|) vs zeta, and the semilocal(p=2) SHIFT, vs truncation K")
for K in [40,80,160,240]:
    a0,b0=jacobi(w0,L,M,K); a2,b2=jacobi(lambda s:wp(s,2),L,M,K)
    e0=np.sort(np.linalg.eigvalsh(Wmat(a0,b0,K)))
    e2=np.sort(np.linalg.eigvalsh(Wmat(a2,b2,K)))
    neg0=e0[e0<0]; neg2=e2[e2<0]
    g0=2*np.sqrt(-neg0[::-1][:5]) if len(neg0)>=5 else np.array([])   # largest |neg| -> ... try both ends
    # gamma from the least-negative (closest to 0) negative eigenvalues:
    g0b=2*np.sqrt(-neg0[-5:][::-1]) if len(neg0)>=5 else np.array([])
    g2b=2*np.sqrt(-neg2[-5:][::-1]) if len(neg2)>=5 else np.array([])
    print(f"\nK={K}: #neg(base)={len(neg0)} #neg(semi)={len(neg2)}")
    if len(g0b)>=3:
        print(f"  base gamma (least-neg 5): {np.round(g0b,3)}")
        print(f"  semi gamma (least-neg 5): {np.round(g2b,3)}")
        m=min(len(g0b),len(g2b))
        print(f"  SHIFT gamma_semi-gamma_base: {np.round(g2b[:m]-g0b[:m],4)}")
    print(f"  base most-neg 5 gamma: {np.round(2*np.sqrt(-neg0[:5]),3) if len(neg0)>=5 else 'NA'}")
