"""If the Weil form Q saturates (near-null space), WHAT are the null modes? Map them to functions on
[-a,a] (sine basis v_k) and characterize. In the Weil/explicit-formula picture the form's degeneracy
should relate to the zeros / the reproducing structure. Also: how does dim(near-null) scale with K,a?"""
import numpy as np, sys
sys.path.insert(0,'operators')
from suzuki_exact_Q import build_Q
for a in [1.0,1.5,2.0]:
    K=max(70,int(30*a)+40)
    Q,_=build_Q(a,K,NG=400,return_pieces=True)
    Q=(Q+Q.T)/2
    w,V=np.linalg.eigh(Q)
    scale=np.median(np.abs(w))
    nnull=int(np.sum(np.abs(w)<1e-6*scale))
    smallpos=int(np.sum((w>1e-6*scale)&(w<1e-2*scale)))
    print(f"a={a} K={K}: eig scale~{scale:.2f}; dim(|λ|<1e-6·scale)={nnull}; #eig in (0, 1e-2·scale)={smallpos}")
    print(f"   smallest 8 eigenvalues: {np.array2string(w[:8], precision=2, floatmode='maxprec')}")
    # rank
    rank=int(np.sum(np.abs(w)>1e-9*scale))
    print(f"   numerical rank = {rank}/{K}  ({'FULL' if rank==K else f'deficient by {K-rank}'})")
    # map the smallest-|λ| eigenvector to a function on [-a,a], sample
    v=V[:,np.argmin(np.abs(w))]
    xs=np.linspace(-a,a,9)
    def vk(k,x): return np.sin((k)*np.pi*(x+a)/(2*a))  # v_k on [-a,a]
    fx=[sum(v[k-1]*vk(k,x) for k in range(1,K+1)) for x in xs]
    print(f"   null-vector-as-function samples on [-a,a]: {np.array2string(np.array(fx),precision=3,floatmode='maxprec')}")
    # dominant sine components of the null vector
    top=np.argsort(-np.abs(v))[:6]+1
    print(f"   null vector dominant sine modes k: {top.tolist()}  weights {np.round(v[top-1],3).tolist()}")
    print()
