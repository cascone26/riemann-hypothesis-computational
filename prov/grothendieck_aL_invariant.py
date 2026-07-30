import numpy as np
from scipy.integrate import solve_ivp
Z=np.load('operators/zeros500.npy')
# rho_op(gamma) = |dphi/dE|/pi  (Riemann's density-of-states of D at the true zero).
# Operator freq knob w: potential FP=w^2, tail phase ph = w*x - E^2/(2 w x).  build w=2pi.
# Fit rho_op(gamma) = a*ln(gamma) + c ; L = exp(-c/a).  TEST: a*L invariant under w-rescale?
#   arithmetic target: rho_true=(1/2pi)ln(gamma/2pi) => a_true=1/2pi=0.15915, L_true=2pi=6.2832,
#                       a_true * L_true = 1.0000  <-- the number that must be hit.
def rho_op(E,X,w):
    mu=-E*E; FP=w*w
    def rhs(x,y): return [y[1],(2*x*y[1]+(FP*x*x-mu)*y[0])/(1-x*x)]
    def ph1(EE):
        s=solve_ivp(rhs,[1+1e-4,X],[1.0,(mu-FP)/2.0] if False else [1.0,(-EE*EE-FP)/2.0],
                    rtol=1e-11,atol=1e-13,max_step=0.01,
                    t_eval=np.linspace(X-5,X,700))
        xs=s.t; u=xs*s.y[0]; p=w*xs-EE*EE/(2*w*xs)
        M=np.column_stack([np.sin(p),np.cos(p)]); (A,B),*_=np.linalg.lstsq(M,u,rcond=None)
        return np.mod(np.arctan2(B,A),np.pi)
    dE=0.01; d=ph1(E+dE)-ph1(E-dE); d=(d+np.pi/2)%np.pi-np.pi/2
    return abs(d/(2*dE))/np.pi
def fit(w,X=170.0,idx=range(0,14)):
    g=np.array([Z[i] for i in idx]); r=np.array([rho_op(gg,X,w) for gg in g])
    a,c=np.polyfit(np.log(g),r,1); L=np.exp(-c/a); return a,c,L,a*L
twopi=2*np.pi
print("  w/2pi     a        c        L=e^-c/a    a*L      [target a*L=1.000]")
for f in [0.5,1.0,2.0]:
    a,c,L,aL=fit(f*twopi)
    print(f"   {f:4.2f}   {a:7.4f}  {c:7.4f}   {L:8.4f}   {aL:7.4f}")
print("\nIf a and L slide with w but a*L holds ~const => a*L is the dilation-invariant obstruction;")
print("no scalar reaches a*L=1. Build a*L vs 1.0 is the archimedean-conductor deficit (my lane).")
