# Higher LP hierarchy: LAGUERRE inequalities L_n[Xi](x)>=0 (pointwise, necessary for RH). Uses the
# Phi-integral: Xi_s(x)=2 int Phi(u)cos(xu)du (scaled, zeros at 2*gamma), Xi^(n)(x)=2 int Phi u^n cos(xu+n*pi/2)du.
# L_1 = f'^2 - f f''.   L_2 = 3 f''^2 - 4 f' f''' + f f''''  (Laguerre, up to positive scaling).
import numpy as np
def Phi(u,nmax=50):
    u=np.atleast_1d(np.asarray(u,float)); n=np.arange(1,nmax+1)[:,None]
    T=(2*np.pi**2*n**4*np.exp(9*u)[None,:]-3*np.pi*n**2*np.exp(5*u)[None,:])*np.exp(-np.pi*n**2*np.exp(4*u)[None,:])
    return T.sum(0)
U=2.2; M=120000; u=np.linspace(1e-7,U,M); du=u[1]-u[0]; phi=Phi(u)
def Xd(x,n):  # n-th derivative of Xi_s at x
    return 2*np.sum(phi*u**n*np.cos(x*u+n*np.pi/2))*du
def L1(x):
    f,f1,f2=Xd(x,0),Xd(x,1),Xd(x,2); return f1*f1-f*f2
def L2(x):
    f,f1,f2,f3,f4=[Xd(x,k) for k in range(5)]; return 3*f2*f2-4*f1*f3+f*f4
xs=np.linspace(0.5,50,400)
l1=np.array([L1(x) for x in xs]); l2=np.array([L2(x) for x in xs])
# normalize by scale to judge tightness: compare min to typical magnitude
print("LAGUERRE inequalities L_n[Xi](x) >= 0 (necessary for RH / LP membership):")
print(f"  L_1: min={l1.min():.4e}  at x={xs[l1.argmin()]:.2f}   (min>=0 => holds; #negative={np.sum(l1<0)})")
print(f"  L_2: min={l2.min():.4e}  at x={xs[l2.argmin()]:.2f}   (min>=0 => holds; #negative={np.sum(l2<0)})")
# tightness: L1 normalized by f'^2 (at zeros L1=f'^2>0); look at ratio L1/(f'^2+f^2) min
def tight1(x):
    f,f1,f2=Xd(x,0),Xd(x,1),Xd(x,2); d=f1*f1+f*f  # scale
    return (f1*f1-f*f2)/d if d>0 else 1.0
t1=np.array([tight1(x) for x in xs])
print(f"  L_1 normalized min (L1/(f'^2+f^2)) = {t1.min():.4f} at x={xs[t1.argmin()]:.2f}  (->0 = knife-edge, >0 comfortable)")
