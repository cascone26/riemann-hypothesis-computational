"""
Step (b) convergence test: is the Weil min eigenvector θ_x captured by the span of the low EVEN prolate functions
{ps_0, ps_2} (Connes' k_λ = h_0,h_4 combo)? Compute residual of θ_x projected onto span{ps_0,ps_2} vs a.
If residual → 0 as a grows: θ_x → k_λ structure → transform → Ξ (step b on track).
"""
import numpy as np
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def weil_min(a, M):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; diffs=np.arange(-M,M+1)*h
    gv=np.array([screw_g(d,Lam,Nmax) for d in diffs]); idx=np.arange(M+1); G=gv[(idx[:,None]-idx[None,:])+M]
    w=np.sqrt(np.full(M+1,h)); w[0]=w[-1]=np.sqrt(0.5*h)
    Gs=w[:,None]*G*w[None,:]; Gs=(Gs+Gs.T)/2
    ev,V=np.linalg.eigh(Gs); e=V[:,0]/w
    return x,h,e/np.sqrt(np.sum(e*e)*h)

def even_prolates(a, x, h, c_bw, k=2):
    # top even prolate spheroidal fns on [-a,a] via sinc concentration kernel, even ones
    D=x[:,None]-x[None,:]
    K=np.where(np.abs(D)<1e-12, c_bw/np.pi, np.sin(c_bw*D)/(np.pi*D))*h
    ev,V=np.linalg.eigh(K)
    evens=[]
    for j in range(len(ev)-1,-1,-1):
        v=V[:,j]
        if np.max(np.abs(v-v[::-1]))<1e-6*np.max(np.abs(v)):   # even
            evens.append(v/np.sqrt(np.sum(v*v)*h))
            if len(evens)>=k: break
    return evens

print("Step(b): residual of Weil θ_x projected onto span{ps_0, ps_2} (low even prolates), vs a")
print(f"{'a':>5} {'1-|<θ,ps0>|²':>13} {'resid after 2 prolates':>22}")
for a in [0.5,1.0,1.5,2.0,3.0]:
    M=min(1000,max(500,int(350/a)))
    x,h,theta=weil_min(a,M)
    # best bandwidth: use c*a~4 (cycle 85); scan a little
    best=None
    for cb in np.linspace(2,12,20):
        ps=even_prolates(a,x,h,cb,2)
        if len(ps)<2: continue
        # gram-schmidt orthonormal (they're ~orthogonal already)
        P=np.array(ps)
        coef=P@ (theta*h)
        proj=coef@P
        resid=np.sqrt(max(0,np.sum((theta-proj)**2)*h))
        ov0=abs(np.sum(theta*ps[0]*h))
        if best is None or resid<best[0]: best=(resid,1-ov0**2)
    print(f"{a:5.2f} {best[1]:13.3e} {best[0]:22.3e}")
print("\nIf residual shrinks with a: θ_x increasingly lives in span of low even prolates (Connes k_λ) → step(b) on track.")
