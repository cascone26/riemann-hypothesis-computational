"""
Exact Schur-complement threshold for kappa=1 + rank-1 lift.
Q_full = Q_0 + P_pole, P_pole rank-2 signature(+,-). Even sector: pole positive part p=sqrt(s+)u+ (cosh);
Q_0 has one negative square phi_1 (even). Schur:  Q_full,even PSD  <=>  s := mu_1 + p_1^2/(1+beta) >= 0,
  p_k=<p,phi_k>,  beta = sum_{k>=2} p_k^2/mu_k.  Equivalently  <p,phi_1>^2/|mu_1| >= 1 + beta.
Verify sign(s) == sign(lambda_a) and that the threshold is ~saturated (knife-edge).
"""
import numpy as np
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def boundary_kernel(t): return -4.0*(np.exp(t/2)+np.exp(-t/2)-2.0)

def analyze(a, M, K):
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
    Qf=(Wp.T@Gf@Wp)/a; Qf=(Qf+Qf.T)/2
    Pp=(Wp.T@Gb@Wp)/a; Pp=(Pp+Pp.T)/2
    Q0=Qf-Pp
    lam_a=np.linalg.eigvalsh(Qf)[0]
    # positive (even) pole direction
    pe,pv=np.linalg.eigh(Pp); s_plus=pe[-1]; u_plus=pv[:,-1]
    p=np.sqrt(max(s_plus,0))*u_plus
    # Q_0 eigendecomp
    mu,phi=np.linalg.eigh(Q0)
    pk=phi.T@p                      # <p,phi_k>
    mu1=mu[0]; p1=pk[0]
    # beta over k>=2 with mu_k>0 (skip near-zero to avoid blowup? report both)
    posmask=mu>1e-9
    beta=np.sum((pk[1:]**2)[mu[1:]>1e-12]/mu[1:][mu[1:]>1e-12])
    s=mu1 + p1*p1/(1+beta)
    LHS=p1*p1/abs(mu1); RHS=1+beta
    return lam_a, mu1, p1*p1, beta, s, LHS, RHS

if __name__=="__main__":
    print("Schur threshold: s=mu_1+p1^2/(1+beta) should match sign(lambda_a); LHS=p1^2/|mu1| vs RHS=1+beta")
    print(f"{'a':>5} {'lambda_a':>11} {'mu_1':>9} {'p1^2':>9} {'beta':>10} {'s(Schur)':>11} {'LHS':>9} {'RHS':>9} {'LHS-RHS':>9}")
    for a in [0.5,1.0,2.0,3.0]:
        M=min(2600,max(1500,int(400/a)))
        la,mu1,p1sq,beta,s,L,R=analyze(a,M=M,K=80)
        print(f"{a:5.2f} {la:11.3e} {mu1:9.4f} {p1sq:9.4f} {beta:10.3f} {s:11.3e} {L:9.4f} {R:9.4f} {L-R:9.3e}",flush=True)
    print("\nIf s~lambda_a (sign+size) and LHS~RHS (saturated): the Schur/BS threshold IS the RH criterion here.")
