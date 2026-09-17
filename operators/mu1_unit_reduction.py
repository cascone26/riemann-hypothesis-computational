import numpy as np, mpmath as mp
mp.mp.dps=25
C0=float(mp.euler)
def kernel(u):  # |u| log|u|, k(0)=0
    au=np.abs(u); out=np.zeros_like(au)
    nz=au>1e-300; out[nz]=au[nz]*np.log(au[nz]); return out
def lambda_min_unit(M,K):
    # deriv form of kernel |t|log|t| on (-1,1), sine-Galerkin, min generalized eig
    x=np.linspace(-1,1,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gv=kernel(diffs); idx=np.arange(M+1); shift=(idx[:,None]-idx[None,:])+M
    G=gv[shift]
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+1),kk)/2.0
    Vp=(kk*np.pi/2.0)*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    Wp=Vp*w[:,None]; Q=Wp.T@G@Wp; Q=(Q+Q.T)/2
    ev=np.linalg.eigvalsh(Q)/1.0   # ||v||^2 = 1*sum c^2 (a=1)
    return ev[0]
print("VERIFY reduction  mu1 = 1 + (1/2)*lambda_min^(1)  [lambda_min^(1) = |t|log|t| deriv form on (-1,1)]",flush=True)
print("a->0 extrapolation gave mu1 ~ 0.2348 (quad) / 0.2434 (linear); predicts lambda_min^(1) ~ -1.530",flush=True)
print(f"{'M':>6} {'K':>4} {'lambda_min^(1)':>15} {'mu1=1+lam/2':>13}",flush=True)
for M,K in [(3000,60),(5000,90),(7000,120),(9000,150)]:
    lm=lambda_min_unit(M,K); mu1=1+0.5*lm
    print(f"{M:6d} {K:4d} {lm:15.6f} {mu1:13.6f}",flush=True)
# Richardson-ish: largest M,K is best
lm=lambda_min_unit(9000,180); mu1=1+0.5*lm
print(f"\nBEST (M=9000,K=180): lambda_min^(1)={lm:.6f}  mu1={mu1:.6f}",flush=True)
# test closed forms for lambda_min^(1)
L2=float(mp.log(2))
for name,val in sorted({"-(1+log2)/... ":0,"-3/2":-1.5,"-(1+C0)":-(1+C0),"-log2-... ":0,
    "-(1+log2)":-1-L2,"-pi^2/... ":0,"-2+log2":-2+L2,"-(2-log2)... ":0,
    "-(1+2log2-... )":0,"C0-2":C0-2,"-1-C0/... ":0,"-2C0-... ":0,
    "-(3-log4)":-(3-2*L2),"-(1+C0+... )":0}.items(),key=lambda kv:abs(kv[1]-lm)):
    if val!=0: print(f"  lam cand {name:16s}={val:.5f} |diff|={abs(val-lm):.5f}",flush=True)
