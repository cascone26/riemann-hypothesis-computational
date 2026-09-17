"""
Connes' EXACT construction: k(u)=E(h)(u)=u^{1/2} sum_{n>=1} h(nu), h(u)=(pi/2)u^2(2pi u^2-3)e^{-pi u^2}.
Riemann/Connes: Fourier(k) = Xi, so zeros of K(t)=int k(u) u^{it} du/u are the zeta zeros.
VERIFY (Handle): K(t) real-axis zeros == zeta zeros. Then this is the exact k_lambda=inf object θ_x approximates.
"""
import numpy as np, mpmath as mp

def h(u): return (np.pi/2)*u**2*(2*np.pi*u**2-3)*np.exp(-np.pi*u**2)

def k(u, N=2000):
    u=np.atleast_1d(u).astype(float); s=np.zeros_like(u)
    for n in range(1,N+1):
        term=h(n*u)
        s+=term
        if np.max(np.abs(term))<1e-18: break
    return np.sqrt(u)*s

# vanishing integral check
from scipy.integrate import quad
I_h,_=quad(lambda u: h(u), 0, 10)
print(f"int_0^inf h du = {I_h:.3e}  (E(h) well-defined; h has ~vanishing full-line integral)")

# K(t) = int_0^inf k(u) u^{it} du/u  = int k(u) u^{it-1} du ; even in log => real for real t
uu=np.linspace(1e-3, 6, 6000); ku=k(uu)
def K(t): return np.trapezoid(ku*np.cos(t*np.log(uu))/uu, uu)   # real part (k even in log-var about u=1? use cos)
ts=np.linspace(1,45,5000); Kv=np.array([K(t) for t in ts])
sign=np.sign(Kv); cross=ts[:-1][np.diff(sign)!=0]
zz=[float(mp.im(mp.zetazero(n))) for n in range(1,11)]
print(f"K(t) zero crossings (first 12): {np.round(cross[:12],4)}")
print(f"zeta zeros (first 10):          {np.round(zz,4)}")
print("match check:")
for g in zz[:8]:
    near=cross[np.argmin(np.abs(cross-g))]
    print(f"   zeta={g:.4f}  K-zero={near:.4f}  diff={abs(near-g):.4f}")
