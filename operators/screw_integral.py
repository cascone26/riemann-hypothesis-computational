"""RH-critical reduction: min-eig(a) ~ (pi^2/4a) integral(g_screw) for the smoothest mode.
Compute integral g_screw(u) du and its sign. If ~0, the next-order term carries RH."""
import numpy as np, os
HERE=os.path.dirname(os.path.abspath(__file__))
src=open(os.path.join(HERE,"suzuki_screw_lambda_a.py")).read()
cut=src.find("if __name__"); exec(src[:cut] if cut>0 else src, globals())
# g_screw over a wide range; integrate; also the weighted integral matching Q(phi0) exactly
for U in (5,8,12,16):
    Nmax=int(np.floor(np.exp(U)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    us=np.linspace(-U,U,20001); du=us[1]-us[0]
    g=np.array([screw_g(u,Lam,Nmax) for u in us])
    I=np.trapezoid(g,us)
    # the actual Q(phi0)-relevant integral: int g(u)(1-|u|/2a)cos(pi u/2a) with a=U/2
    a=U/2.0
    wt=np.where(np.abs(us)<2*a,(1-np.abs(us)/(2*a))*np.cos(np.pi*us/(2*a)),0)
    Iw=np.trapezoid(g*wt,us)
    print(f"U={U}: int g_screw = {I:+.5f} | Q-weighted int(a={a}) = {Iw:+.5f} | g(0)={g[len(g)//2]:.4f}", flush=True)
print("DONE")
