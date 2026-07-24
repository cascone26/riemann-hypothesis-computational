# Krein-Langer nonneg-cosine CERTIFICATE for lambda_a>=0. Fit g(t) on [0,2a] by  -c2*t^2 + c + sum m_j cos(xi_j t),
# c2,m_j>=0 (NNLS). Exact fit => Q(v)=sum m_j xi_j^2 |vhat(xi_j)|^2 >=0 for all v in (-a,a) => lambda_a>=0 PROVEN.
# Watch: residual->0? atoms -> true zeros gamma=14.13,21.02,25.01,30.42,32.94?
import numpy as np, sys, os
from scipy.optimize import nnls
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
import suzuki_screw_lambda_a as S
G_TRUE=[14.134725,21.022040,25.010858,30.424876,32.935062,37.586178,40.918719,43.327073]
for a in [0.37, 0.48, 0.60, 0.80, 1.0]:
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
    ts=np.linspace(0,2*a,2000)
    g=np.array([S.screw_g(t,Lam,Nmax) for t in ts])
    # atom grid: dense near expected zeros + broad coverage
    xis=np.unique(np.concatenate([np.linspace(0.2,60,3000), np.linspace(60,500,1500),
                                  np.concatenate([np.linspace(gg-0.4,gg+0.4,40) for gg in G_TRUE])]))
    cols=[ -ts**2, np.ones_like(ts), -np.ones_like(ts) ] + [np.cos(x*ts) for x in xis]
    A=np.column_stack(cols)
    m,res=nnls(A,g,maxiter=60000)
    relres=res/np.linalg.norm(g)
    ma=m[3:]; nat=int(np.sum(ma>1e-9))
    ordr=np.argsort(ma)[::-1][:6]; tops=[(round(float(xis[i]),2),round(float(ma[i]),4)) for i in ordr if ma[i]>1e-9]
    print(f"a={a:.2f}: rel-res={relres:.3e}  c2={m[0]:.4f}  atoms={nat}  top(xi,m)={tops}", flush=True)
print("\nTrue zeros: 14.13,21.02,25.01,30.42,32.94 | residual->0 & atoms->zeros => Krein certificate viable for first step.")
