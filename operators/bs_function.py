"""
Birman-Schwinger verification: Q_full = Q_0 + P_pole (P_pole = pole/boundary Galerkin matrix).
(1) rank of P_pole (expect 2: e^{x/2},e^{-x/2}; 1 on even sector).
(2) BS(lam) = generalized: lowest eig of Q_full should satisfy the rank-one Krein relation.
    We verify via the exact identity: eigenvalues of Q_full where P_pole=sum_j s_j |p_j><p_j| (s_j>0).
    Check: det(Q_0 - lam) * det(I + [<p_i|(Q_0-lam)^{-1}|p_j> s_j]) = 0 at lam=lambda_a.
"""
import numpy as np
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def boundary_kernel(t): return -4.0*(np.exp(t/2)+np.exp(-t/2)-2.0)

def build(a, M, K):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gf=np.array([screw_g(d,Lam,Nmax) for d in diffs]); gb=boundary_kernel(diffs)
    idx=np.arange(M+1); shift=(idx[:,None]-idx[None,:])+M
    Gf=gf[shift]; Gb=gb[shift]
    kk=np.arange(1,K+1)
    arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    Wp=Vp*w[:,None]
    Qf=(Wp.T@Gf@Wp); Qf=(Qf+Qf.T)/2
    Ppole=(Wp.T@Gb@Wp); Ppole=(Ppole+Ppole.T)/2
    Q0=Qf-Ppole
    return Qf/a, Q0/a, Ppole/a

if __name__=="__main__":
    for a in [1.0, 2.0]:
        M=min(2600,max(1500,int(400/a)))
        Qf,Q0,P=build(a,M=M,K=80)
        # rank of pole matrix
        pe=np.linalg.eigvalsh(P)
        rankP=np.sum(np.abs(pe)>1e-6*max(1,abs(pe).max()))
        npos=np.sum(pe>1e-6*max(1,abs(pe).max())); nneg=np.sum(pe<-1e-6*max(1,abs(pe).max()))
        lam_full=np.linalg.eigvalsh(Qf)[0]
        e0=np.linalg.eigvalsh(Q0)
        mu1=e0[0]
        # BS relation: at lam=lam_full, Q0 - lam_full is invertible (mu1<0<... lam~0>mu1),
        # and det(I + (Q0-lam)^{-1} P) = 0. Check smallest singular value of (I + (Q0-lam)^{-1}P).
        lam=lam_full
        R=np.linalg.solve(Q0-lam*np.eye(len(Q0)), P)   # (Q0-lam)^{-1} P
        Mmat=np.eye(len(Q0))+R
        sv=np.linalg.svd(Mmat, compute_uv=False)
        print(f"a={a}: rank(P_pole)={rankP} (pos {npos}, neg {nneg}); lam_a={lam_full:.3e}, mu1(Q0)={mu1:.4f}")
        print(f"   BS check: min singular value of (I+(Q0-lam_a)^-1 P) = {sv.min():.3e}  (=0 confirms rank-perturb reln)")
        # RH-threshold form: does the ONE negative eigenvector phi1 of Q0 carry the pole overlap?
        w1=np.linalg.eigh(Q0)[1][:,0]
        # pole vector(s): P = sum s_j p_j p_j^T ; dominant pole direction
        pev,pvec=np.linalg.eigh(P)
        pdir=pvec[:,-1]  # largest pole eigenvector
        ov=abs(np.dot(w1,pdir))
        print(f"   |<phi1(neg eigvec), dominant pole dir>| = {ov:.4f}  (pole must overlap the neg square to lift it)")
