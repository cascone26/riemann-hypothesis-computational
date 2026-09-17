import numpy as np, sys, os
sys.path.insert(0,'/Users/scones/projects/riemann/operators'); os.chdir('/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto
def boundary_kernel(t): return -4.0*(np.exp(t/2)+np.exp(-t/2)-2.0)
def build(a,M,K):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; diffs=np.arange(-M,M+1)*h
    gf=np.array([screw_g(d,Lam,Nmax) for d in diffs]); gb=boundary_kernel(diffs)
    idx=np.arange(M+1); shift=(idx[:,None]-idx[None,:])+M
    Gf=gf[shift]; Gb=gb[shift]; kk=np.arange(1,K+1)
    arg=np.pi*np.outer((x+a),kk)/(2*a); Vp=(kk*np.pi/(2*a))*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h; Wp=Vp*w[:,None]
    Qf=(Wp.T@Gf@Wp)/a; Qf=(Qf+Qf.T)/2; Pp=(Wp.T@Gb@Wp)/a; Pp=(Pp+Pp.T)/2
    return Qf, Qf-Pp, Pp
print("a     M(LHS-RHS)   -(1+<p,Q0^-1 p>)   min-eig(Q_full)   sign-match")
for a in [0.5,1.0,2.0,3.0]:
    Mgrid=min(2600,max(1500,int(400/a)))
    Qf,Q0,Pp=build(a,Mgrid,80)
    pe,pv=np.linalg.eigh(Pp); p=np.sqrt(max(pe[-1],0))*pv[:,-1]
    mu,phi=np.linalg.eigh(Q0); pk=phi.T@p; mu1=mu[0]; p1=pk[0]
    beta=np.sum((pk[1:]**2)[mu[1:]>1e-12]/mu[1:][mu[1:]>1e-12])
    LHS=p1*p1/abs(mu1); RHS=1+beta; M=LHS-RHS
    pQ0inv_p=float(p@np.linalg.solve(Q0,p))     # <p, Q0^-1 p>
    ident=-(1+pQ0inv_p)
    mineig_full=np.linalg.eigvalsh(Qf).min()
    ok = (np.sign(M)==np.sign(ident)) and abs(M-ident)/max(abs(M),1e-30)<0.05
    print(f"{a:4.2f}  {M:11.4e}  {ident:15.4e}   {mineig_full:14.4e}   {'YES' if ok else 'NO'}",flush=True)
print("Identity M = -(1+<p,Q0^-1 p>) confirms: RH-corner <=> rank-1 pole lifts Q0's single neg eigenvalue above 0.")
