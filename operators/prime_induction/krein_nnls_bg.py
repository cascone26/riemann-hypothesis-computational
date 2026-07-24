import numpy as np, sys, os, time
from scipy.optimize import nnls
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
import suzuki_screw_lambda_a as S
G_TRUE=[14.134725,21.022040,25.010858,30.424876,32.935062,37.586178,40.918719,43.327073,48.005151,49.773832]
def fit(a, natom_broad=1200):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
    ts=np.linspace(0,2*a,1400); g=np.array([S.screw_g(t,Lam,Nmax) for t in ts])
    xis=np.unique(np.concatenate([np.linspace(0.2,80,natom_broad), np.linspace(80,260,300),
                    np.concatenate([np.linspace(gg-0.3,gg+0.3,25) for gg in G_TRUE])]))
    A=np.column_stack([-ts**2, np.ones_like(ts), -np.ones_like(ts)]+[np.cos(x*ts) for x in xis])
    m,res=nnls(A,g,maxiter=40000); relres=res/np.linalg.norm(g)
    ma=m[3:]; ordr=np.argsort(ma)[::-1][:6]
    tops=[(round(float(xis[i]),2),round(float(ma[i]),4)) for i in ordr if ma[i]>1e-9]
    return relres, int(np.sum(ma>1e-9)), float(m[0]), tops
with open('/tmp/pi/krein_out.txt','w') as f:
    for a in [0.37,0.48,0.60,0.80,1.0]:
        t0=time.time(); rr,nat,c2,tops=fit(a)
        line=f"a={a:.2f}: rel-res={rr:.3e}  c2={c2:.4f}  atoms={nat}  top={tops}  ({time.time()-t0:.0f}s)"
        f.write(line+"\n"); f.flush()
    f.write("DONE\n")
