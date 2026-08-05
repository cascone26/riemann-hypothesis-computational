"""Clean small-a mu1 extraction (RELIABLE regime a<0.347, purely archimedean). Extract mu1 to
high precision from lambda_a(a) = c log(1/a) + b + e a + ..., mu1 = b - (log2pi - psi(2))/2-ish.
Then hunt closed form."""
import numpy as np, mpmath as mp, os
HERE=os.path.dirname(os.path.abspath(__file__))
src=open(os.path.join(HERE,"suzuki_screw_lambda_a.py")).read()
cut=src.find("if __name__"); exec(src[:cut] if cut>0 else src, globals())
avals=np.array([0.05,0.08,0.11,0.14,0.18,0.22,0.26,0.30])
las=[]
for a in avals:
    la=lambda_a(a, M=8000, K=60)
    las.append(la)
    print(f"a={a:.2f}: lambda_a={la:.6f}", flush=True)
las=np.array(las)
# fit c log(1/a) + b + e a
X=np.column_stack([np.log(1/avals), np.ones_like(avals), avals])
coef,*_=np.linalg.lstsq(X,las,rcond=None)
c,b,e=coef
const=float(0.5*(mp.log(2*mp.pi)-mp.digamma(2)))
mu1=b-const
print(f"\nfit: c={c:.4f} (Suzuki Thm1.4 =1), b={b:.4f}, e={e:.3f}")
print(f"const=(log2pi-psi(2))/2={const:.4f}; mu1=b-const={mu1:.4f}")
# closed-form candidates near mu1
print("closed-form candidates:")
for name,val in [("gamma_E/2",mp.euler/2),("1/4",0.25),("(1-gamma)/2",(1-mp.euler)/2),
                 ("gamma_E",mp.euler),("log2/e",mp.log(2)/mp.e),("1-gamma",1-mp.euler),
                 ("pi^2/40",mp.pi**2/40),("gamma^2",mp.euler**2)]:
    print(f"  {name} = {float(val):.4f}  diff={abs(float(val)-mu1):.4f}")
print("DONE")
