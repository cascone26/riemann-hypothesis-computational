# KREIN-LANGER CERTIFICATE via NNLS: fit g(t) on [0,2a] by  c + sum_j m_j cos(xi_j t),  m_j>=0, c free.
# Exact fit => Q(v) = sum m_j xi_j^2 |vhat(xi_j)|^2 >= 0 for ALL v supported in (-a,a)  => lambda_a >= 0 PROVEN.
# (Under RH the global version has atoms at the zeros, m=2/gamma^2 — §47. This is the finite-window analog.)
import numpy as np, sys, os
from scipy.optimize import nnls
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
import suzuki_screw_lambda_a as S
for a in [0.30, 0.37, 0.42, 0.48, 0.54, 0.60, 0.80]:
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
    ts=np.linspace(0,2*a,1200)
    g=np.array([S.screw_g(t,Lam,Nmax) for t in ts])
    xis=np.concatenate([np.linspace(0.3,60,900), np.linspace(60.2,300,600)])
    # design: [1, cos(xi_j t)]; c free => allow +/- via two columns
    A=np.column_stack([np.ones_like(ts), -np.ones_like(ts)] + [np.cos(x*ts) for x in xis])
    m,res=nnls(A,g,maxiter=10000)
    relres=res/np.linalg.norm(g)
    nat=np.sum(m[2:]>1e-9)
    # top atoms
    ord_=np.argsort(m[2:])[::-1][:5]; tops=[(round(float(xis[i]),2),round(float(m[2:][i]),4)) for i in ord_ if m[2:][i]>1e-9]
    print(f"a={a:.2f}: NNLS rel-residual={relres:.2e}  atoms={nat}  top(xi,m)={tops}", flush=True)
print("\nrel-residual ~ 0 => nonneg-cosine certificate EXISTS on [0,2a] => lambda_a>=0 proven at that a.")
print("watch: do top atoms approach the true zeros gamma=14.13,21.02,... as a grows? (spectral convergence)")
