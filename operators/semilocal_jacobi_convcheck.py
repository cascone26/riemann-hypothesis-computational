import numpy as np
from scipy.special import loggamma
def w0(s): return np.exp(2.0*np.real(loggamma(0.25+0.5j*s)))
def wp(s,p): return w0(s)/(1.0-2.0*p**-0.5*np.cos(s*np.log(p))+1.0/p)
def jac(weight,L,M,K=14):
    s=np.linspace(-L,L,M); ds=s[1]-s[0]; w=weight(s)*ds
    Mt=w.sum(); qp=np.zeros_like(s); q=np.sqrt(w/Mt); beta=np.zeros(K); b=0.0
    for k in range(K):
        v=s*q; a=np.dot(q,v); v=v-a*q-b*qp; v=v-q*np.dot(q,v)-qp*np.dot(qp,v)
        b=np.linalg.norm(v)
        if k+1<K: beta[k+1]=b
        qp=q; q=v/b
    return beta
print("Convergence of delta b_n (p=2) across (L,M) grids -- is the imprint grid-stable?")
cfgs=[(35,120000),(45,200000),(55,300000)]
tabs={}
for L,M in cfgs:
    tabs[(L,M)]=jac(lambda s:wp(s,2),L,M)-jac(w0,L,M)
hdr="  n"+"".join(f"   L={L},M={M//1000}k" for L,M in cfgs)
print(hdr)
for n in range(1,12):
    print(f"{n:3d}"+"".join(f"{tabs[c][n]:15.5f}" for c in cfgs))
spread=max(abs(tabs[cfgs[0]][n]-tabs[cfgs[2]][n]) for n in range(1,10))
print(f"\nmax|db(L=35) - db(L=55)| over n=1..9 = {spread:.2e}  (converged if << imprint size ~1)")
