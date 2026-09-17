import numpy as np
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto
def boundary_kernel(t): return -4.0*(np.exp(t/2)+np.exp(-t/2)-2.0)
def build_Q(a,M,K):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gf=np.array([screw_g(d,Lam,Nmax) for d in diffs])
    idx=np.arange(M+1); shift=(idx[:,None]-idx[None,:])+M
    G=gf[shift]
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    Wp=Vp*w[:,None]; Q=Wp.T@G@Wp; return (Q+Q.T)/2
print("AF positive-index/rank-trace lens on OUR localized Weil (screw) form Q_full:")
print(f"{'a':>5} {'K':>4} {'posIdx':>7} {'negIdx':>7} {'trQ':>10} {'||Q||_HS^2':>12} {'(trP+)^2/||P+||^2':>16} {'posIdx (CS bnd)':>14}")
for a in [0.5,1.0,2.0,3.0]:
    K=80; ev=np.linalg.eigvalsh(build_Q(a,3000,K))
    tol=1e-9*max(1,abs(ev).max())
    pos=ev[ev>tol]; neg=ev[ev<-tol]
    trQ=ev.sum(); hs=(ev**2).sum()
    trPp=pos.sum(); hsPp=(pos**2).sum()
    cs=(trPp**2/hsPp) if hsPp>0 else 0
    print(f"{a:5.2f} {K:4d} {len(pos):7d} {len(neg):7d} {trQ:10.3f} {hs:12.3f} {cs:16.3f} {len(pos):14d}")
print("\nNote: (trP+)^2/||P+||^2 <= posIdx (Cauchy-Schwarz) is the AF rank-lower-bound proxy.")
print("If this ratio tracks a fixed fraction of K, that's an AF-style proportion signal; else objects differ.")
