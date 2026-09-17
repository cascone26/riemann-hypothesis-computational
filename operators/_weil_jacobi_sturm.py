import numpy as np, math
import operators.connes_completeness_highT as C
# CM-DEFERRED OBJECT: the semilocal Jacobi matrix of the cyclic pair (Weil form → orthogonal-poly recurrence).
# Lanczos-tridiagonalize the Weil form tau; get Jacobi coeffs (a_n diag, b_n off-diag). Then Sturm:
# #neg eigenvalues = #sign-changes in leading principal minors D_k. Attack kappa<=1 structurally.

def build_weil(lsq, N):
    L=float(np.log(lsq)); VM=C.vm_table(lsq)
    tau,idx=C.build_tau(N,L,VM)
    return tau

def lanczos(A, v0, m):
    n=A.shape[0]; V=np.zeros((n,m)); a=np.zeros(m); b=np.zeros(m)
    v=v0/np.linalg.norm(v0); vprev=np.zeros(n); bprev=0.0
    for j in range(m):
        V[:,j]=v; w=A@v; a[j]=v@w; w=w-a[j]*v-bprev*vprev
        # full reorth
        w=w-V[:,:j+1]@(V[:,:j+1].T@w)
        bprev=np.linalg.norm(w)
        if bprev<1e-12: b[j]=0; break
        b[j]=bprev; vprev=v; v=w/bprev
    return a[:j+1], b[:j] if j<m else b[:m-1], j+1

def sturm_neg_count(a, b):
    # leading principal minors of Jacobi J: D_{-1}=1,D_0=a_0, D_k=a_k D_{k-1} - b_{k-1}^2 D_{k-2}
    # #neg eigenvalues = #sign changes in sequence 1,D_0,D_1,...  (for eigenvalue 0 threshold)
    Dm2=1.0; Dm1=a[0]; signs=[1.0, np.sign(Dm1) if Dm1!=0 else 1.0]
    for k in range(1,len(a)):
        Dk=a[k]*Dm1 - (b[k-1]**2)*Dm2
        signs.append(np.sign(Dk) if Dk!=0 else signs[-1])
        Dm2,Dm1=Dm1,Dk
    changes=sum(1 for i in range(1,len(signs)) if signs[i]*signs[i-1]<0)
    return changes

for lsq,N in [(14,40),(60,50),(250,60)]:
    tau=build_weil(lsq,N)
    D=tau.shape[0]
    # cyclic vector: the "Gaussian"/uniform vector (natural cyclic vector for the pair)
    v0=np.ones(D)
    a,b,m=lanczos(tau, v0, D)
    nc=sturm_neg_count(a,b)
    ev=np.linalg.eigvalsh(tau); kap_direct=int((ev<-1e-9).sum())
    print(f"lsq={lsq} N={N} dim={D}: Lanczos steps={m}")
    print(f"   Jacobi a_n (diag) first6: {np.round(a[:6],4).tolist()}")
    print(f"   Jacobi b_n (offdiag) first6: {np.round(b[:6],4).tolist()}  (all>0 => genuine Jacobi/positive measure)")
    print(f"   Sturm #neg (sign-changes in minors) = {nc}   direct kappa = {kap_direct}", flush=True)
