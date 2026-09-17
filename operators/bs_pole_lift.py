"""
Birman-Schwinger structural test: is the Weil form WITHOUT the pole/boundary term negative, with the
rank-one pole term lifting it to lambda_a ~ 0?
  g_full = boundary + arch_core + arch_tail + prime.  boundary = -4(e^{t/2}+e^{-t/2}-2) = RANK-ONE (cycle 86).
  Q_0 = int int (g_full - boundary) w'w'  = form without poles.  mu_0(a) = inf spec Q_0.
Compare mu_0(a) to lambda_a(full). If mu_0<0 and full ~0, the pole (rank-one, +) lifts a negative core to ~0.
"""
import numpy as np, time
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def boundary_kernel(t):
    return -4.0*(np.exp(t/2)+np.exp(-t/2)-2.0)

def specs(a, M, K):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gf=np.array([screw_g(d,Lam,Nmax) for d in diffs])
    gb=boundary_kernel(diffs)
    idx=np.arange(M+1); shift=(idx[:,None]-idx[None,:])+M
    Gf=gf[shift]; G0=(gf-gb)[shift]     # full minus boundary
    kk=np.arange(1,K+1)
    arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    Wp=Vp*w[:,None]
    Qf=Wp.T@Gf@Wp; Qf=(Qf+Qf.T)/2
    Q0=Wp.T@G0@Wp; Q0=(Q0+Q0.T)/2
    lam_full=(np.linalg.eigvalsh(Qf)/a)[0]
    mu0=(np.linalg.eigvalsh(Q0)/a)[0]
    return lam_full, mu0

if __name__=="__main__":
    print("BS test: mu_0 = inf spec (Weil form WITHOUT pole/boundary) vs lambda_a (full)")
    print(f"{'a':>5} {'lambda_a(full)':>15} {'mu_0(no poles)':>15} {'pole lift':>12}")
    for a in [0.5,0.75,1.0,1.5,2.0,3.0]:
        M=min(3000,max(1500,int(400/a)))
        lf,m0=specs(a,M=M,K=80)
        print(f"{a:5.2f} {lf:15.6e} {m0:15.6e} {lf-m0:12.5f}",flush=True)
    print("\nif mu_0<0 and lambda_a~0: the rank-one pole term (poles s=0,1) is ESSENTIAL — it lifts the")
    print("negative arch+prime core up to the knife-edge. RH <=> pole lift exactly compensates the core.")
