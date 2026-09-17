import numpy as np
import operators.connes_completeness_highT as C

def build_weil(lsq,N):
    L=float(np.log(lsq)); VM=C.vm_table(lsq); tau,idx=C.build_tau(N,L,VM); return tau

def lanczos_full(A, v0):
    n=A.shape[0]; V=np.zeros((n,n)); a=[]; b=[]
    v=v0/np.linalg.norm(v0); vp=np.zeros(n); bp=0.0
    for j in range(n):
        V[:,j]=v; w=A@v; aj=v@w; a.append(aj)
        w=w-aj*v-bp*vp
        w=w-V[:,:j+1]@(V[:,:j+1].T@w)   # full reorth
        w=w-V[:,:j+1]@(V[:,:j+1].T@w)   # twice
        bj=np.linalg.norm(w)
        if bj<1e-10: break
        b.append(bj); vp=v; v=w/bj
    return np.array(a), np.array(b)

np.random.seed(1)
for lsq,N in [(14,60),(250,80)]:
    tau=build_weil(lsq,N); D=tau.shape[0]
    v0=np.random.randn(D)          # random => cyclic w.p.1
    a,b=lanczos_full(tau,v0)
    # faithfulness: eigenvalues of Jacobi J vs tau
    J=np.diag(a)+np.diag(b,1)+np.diag(b,-1)
    evJ=np.sort(np.linalg.eigvalsh(J)); evT=np.sort(np.linalg.eigvalsh(tau))
    match=np.max(np.abs(evJ-evT[:len(evJ)])) if len(evJ)==len(evT) else 9.9
    kap=int((evT<-1e-9).sum())
    print(f"lsq={lsq} N={N} dim={D}: Lanczos len={len(a)}, |eig(J)-eig(tau)|={match:.2e} (faithful if ~0), kappa={kap}")
    # Remling: do a_n, b_n converge to limits? show tails
    print(f"  a_n tail (last 8): {np.round(a[-8:],3).tolist()}")
    print(f"  b_n tail (last 8): {np.round(b[-8:],3).tolist()}")
    # essential spectrum edge estimate a* - 2 b* from tail means
    astar=np.mean(a[len(a)//2:]); bstar=np.mean(b[len(b)//2:])
    print(f"  tail-mean a*={astar:.3f} b*={bstar:.3f}  ess.spec edge a*-2b*={astar-2*bstar:.3f} (>0 => positive continuum, finite bound states)", flush=True)
