"""
Krein-Langer test: count NEGATIVE eigenvalues of the no-pole Weil form Q_0 vs the full form Q_full.
If Q_0 has exactly k negative squares and the rank-<=2 pole term (e^{x/2},e^{-x/2}) removes them so
Q_full is PSD (0 negatives), that IS the Krein-Langer "kappa negative squares = the poles" structure.
Also report the SECOND eigenvalue of Q_0 (interlacing bound: rank-2 lift can raise lowest to <= 3rd eig).
"""
import numpy as np
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def boundary_kernel(t): return -4.0*(np.exp(t/2)+np.exp(-t/2)-2.0)

def counts(a, M, K):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gf=np.array([screw_g(d,Lam,Nmax) for d in diffs]); gb=boundary_kernel(diffs)
    idx=np.arange(M+1); shift=(idx[:,None]-idx[None,:])+M
    Gf=gf[shift]; G0=(gf-gb)[shift]
    kk=np.arange(1,K+1)
    arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    Wp=Vp*w[:,None]
    Qf=Wp.T@Gf@Wp; Qf=(Qf+Qf.T)/2
    Q0=Wp.T@G0@Wp; Q0=(Q0+Q0.T)/2
    ef=np.linalg.eigvalsh(Qf)/a; e0=np.linalg.eigvalsh(Q0)/a
    tol=1e-4*max(1,abs(e0).max())
    return (ef<-tol).sum(),(e0<-tol).sum(), e0[:4], ef[0]

if __name__=="__main__":
    print("Krein-Langer neg-square count: full Q vs no-pole Q_0")
    print(f"{'a':>5} {'#neg full':>10} {'#neg Q_0':>9} {'Q_0 lowest 4 eigs':>44} {'lam_a':>11}")
    for a in [0.5,1.0,2.0,3.0]:
        M=min(3000,max(1500,int(400/a)))
        nf,n0,e0four,lf=counts(a,M=M,K=80)
        s="["+", ".join(f"{v:.3f}" for v in e0four)+"]"
        print(f"{a:5.2f} {nf:10d} {n0:9d} {s:>44} {lf:11.2e}",flush=True)
    print("\nInterpretation: #neg(Q_0) = kappa = negative squares removed by the rank-<=2 pole term.")
