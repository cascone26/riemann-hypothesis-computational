"""Why EXACTLY κ=2 negative squares appear under prime-damping? Identify the 2 modes that go negative:
are they a symmetric ±pair, pole-related, or tied to specific low zeros? Track the 2 lowest eigenvectors
of Q(δ) as δ:0->small and see what functions they are."""
import numpy as np, sys
sys.path.insert(0,'operators')
from heated_weil import build_Q_heated
a=1.5; K=85
Q0=build_Q_heated(a,K,0.0)
w0,V0=np.linalg.eigh(Q0)
Qd=build_Q_heated(a,K,0.01)
wd,Vd=np.linalg.eigh(Qd)
print(f"a={a}: at δ=0.01, lowest 4 eigenvalues: {np.round(wd[:4]/a,4)}")
print(f"        the 2 negative ones: {np.round(wd[:2]/a,4)}")
# express the 2 negative eigenvectors as functions; check parity (even/odd about x=0)
def vk(k,x,a): return np.sin(k*np.pi*(x+a)/(2*a))
xs=np.linspace(-a,a,11)
for m in range(2):
    v=Vd[:,m]
    fx=np.array([sum(v[k-1]*vk(k,x,a) for k in range(1,K+1)) for x in xs])
    # parity: compare f(x) vs f(-x)
    fxr=fx[::-1]
    par = 'EVEN' if np.allclose(fx,fxr,atol=0.1) else ('ODD' if np.allclose(fx,-fxr,atol=0.1) else 'mixed')
    top=np.argsort(-np.abs(v))[:5]+1
    print(f"   neg mode {m}: parity={par}  dominant sine k={top.tolist()}  (even k=even func? odd k=odd)")
    print(f"      samples: {np.round(fx,2)}")
# overlap of the 2 neg modes with the pole directions E_half_p, E_half_m
from suzuki_exact_Q import E_p
kk=np.arange(1,K+1)
Ehp=E_p(0.5,a,kk); Ehm=E_p(-0.5,a,kk)
Ehp/=np.linalg.norm(Ehp); Ehm/=np.linalg.norm(Ehm)
for m in range(2):
    v=Vd[:,m]
    print(f"   neg mode {m}: |<mode|E_+>|={abs(v@Ehp):.3f}  |<mode|E_->|={abs(v@Ehm):.3f}  (pole overlap)")
