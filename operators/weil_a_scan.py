"""THE RH TEST: converged min-eig(a) vs interval length a. RH <=> stays >0 for all a.
Fit the a-dependence; is it C/a^2 with C>0 (from above, RH-safe) or crossing 0?"""
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
    Q=Wp.T@G@Wp; return (Q+Q.T)/2/a
res=[]
for a in (1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0):
    K=max(40,int(18*a)); M=max(4000,int(1200*a))
    ev=np.linalg.eigvalsh(Qmat(a,M,K))
    res.append((a,ev[0]))
    print(f"a={a:.1f}: min-eig={ev[0]:+.5e} [{'POS' if ev[0]>0 else 'NEG'}]", flush=True)
aa=np.array([r[0] for r in res]); mm=np.array([r[1] for r in res])
print(f"\nall min-eig > 0: {all(mm>0)}")
# fit log(min-eig) ~ log a to get power
if all(mm>0):
    p=np.polyfit(np.log(aa),np.log(mm),1)
    print(f"min-eig ~ a^{p[0]:.2f} (if ~ -2 and stays>0 -> RH-safe saturation from above)")
print("DONE")
