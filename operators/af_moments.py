import sys,os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import ctx
import numpy as np
from dh_certified_tau import build_arb_tau
# AF: proportion >= 2 - R, R = ||Gt||_HS^2 / tr(Gt) = tr(tau^2)/tr(tau) (shared nonzero spectrum).
# 4th-moment Christoffel(m=2): prop >= 1 - Lambda_2(0), with moments (m1,m2,m4). RMT/HL target m_k(1)=1,4/3,2,13/4.
for c in [50,100]:
    ctx.prec=800
    A,DIM=build_arb_tau(c,64,800)
    M=np.array([[float(A[i,j].mid()) for j in range(DIM)] for i in range(DIM)])
    ev=np.linalg.eigvalsh(0.5*(M+M.T)); ev=ev[ev>1e-12]  # positive spectrum (drop razor-nulls)
    d=len(ev)
    m1=ev.sum(); m2=(ev**2).sum(); m3=(ev**3).sum(); m4=(ev**4).sum()
    R=m2/m1
    print(f"c={c}: #pos-eig d={d}, normalized moments tr(t^k)/d: "
          f"m1={m1/d:.4f} m2={m2/d:.4f} m3={m3/d:.4f} m4={m4/d:.4f}")
    print(f"       R = tr(t^2)/tr(t) = {R:.4f}   (AF indicator target 4/3={4/3:.4f})   => AF proportion 2-R = {2-R:.4f}")
    # normalized-to-mean moments (rescale so m1/d=1): mu_k = (tr t^k/d)/(tr t/d)^k
    mu1,mu2,mu4 = 1.0, (m2/d)/(m1/d)**2, (m4/d)/(m1/d)**4
    print(f"       mean-normalized mu2={mu2:.4f} (target 4/3), mu4={mu4:.4f} (HL target 13/4={13/4:.4f})")
print("If mu2~4/3 & mu4~13/4, our finite unconditional instrument REALIZES AF's moments -> path to 13/18 via 4th moment.")
