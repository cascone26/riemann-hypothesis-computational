"""Extract the h^2 coefficient C of min-eig ~ C h^2 (h~1/K). C>0 => saturation from above => RH.
Then analyze the near-null eigenvector: is its Weil value a manifestly positive quantity?"""
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
a=3.0
Ks=np.array([30,40,50,60,70,80,90])
mins=[]
for K in Ks:
    Q,x,kk=Qmat(a,4500,K)
    ev=np.linalg.eigvalsh(Q)
    mins.append(ev[0])
    print(f"K={K}: min-eig={ev[0]:+.5e}", flush=True)
mins=np.array(mins)
# fit min-eig ~ C/K^2 (h = 2a/K): if all >0 and ~1/K^2 with C>0, saturation from above
h=2*a/Ks
b,c=np.polyfit(h**2, mins, 1)
print(f"\nfit min-eig = {c:+.3e} + {b:+.3e} h^2")
print(f"all min-eig > 0: {all(mins>0)}; slope(h^2)>0: {b>0}; intercept(h->0): {c:+.3e}")
print("RH-consistent if all>0 AND intercept ~0 from above")
print("DONE")
