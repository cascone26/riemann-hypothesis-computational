"""RH attack via Weil positivity. Q_W PSD <=> RH. Screw kernel (561-562 machinery).
Push: does min-eig stay STRICTLY >0 as resolution increases (RH holds), and WHY?
Characterize the near-null mode + the positivity mechanism."""
import numpy as np, os
HERE=os.path.dirname(os.path.abspath(__file__))
src=open(os.path.join(HERE,"suzuki_screw_lambda_a.py")).read()
cut=src.find("if __name__"); exec(src[:cut] if cut>0 else src, globals())
def Qmat(a,M,K):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gvals=np.array([screw_g(d,Lam,Nmax) for d in diffs])
    idx=np.arange(M+1); G=gvals[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h; Wp=Vp*w[:,None]
    Q=Wp.T@G@Wp; return (Q+Q.T)/2/a, x, kk
# fixed a, INCREASE K (resolution) -> does min-eig -> 0+ (positive) or cross 0?
a=3.0
print(f"a={a}: min-eig vs resolution K (RH <=> stays >=0):")
prev=None
for K in (20,30,40,50,60,70,80):
    Q,x,kk=Qmat(a,4000,K)
    ev=np.linalg.eigvalsh(Q)
    sign="POS" if ev[0]>0 else "NEG"
    print(f"  K={K:>3}: min-eig={ev[0]:+.5e} [{sign}] 2nd={ev[1]:+.4e} min/max={ev[0]/ev[-1]:+.3e}", flush=True)
print("DONE")
