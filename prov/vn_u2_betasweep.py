import numpy as np
from scipy.integrate import solve_ivp
Z=np.load('operators/zeros500.npy'); FP=(2*np.pi)**2; w=2*np.pi
eps=1e-4; X=120.0
# Two independent solutions at x=1+eps: regular u_r (IC [1,(mu-FP)/2]) and singular u_l (~log(x-1)).
# General x=1 self-adjoint extension: IC = cos(beta)*ICr + sin(beta)*ICl.  beta=0 => Friedrichs (build).
def phase(E,beta):
    mu=-E*E
    def rhs(x,y): return [y[1],(2*x*y[1]+(FP*x*x-mu)*y[0])/(1-x*x)]
    ICr=np.array([1.0,(mu-FP)/2.0])
    ICl=np.array([np.log(eps), 1.0/eps])          # leading singular solution ~ log(x-1)
    IC = np.cos(beta)*ICr + np.sin(beta)*ICl
    xs=np.linspace(X-5,X,600)
    s=solve_ivp(rhs,[1+eps,X],IC,rtol=1e-10,atol=1e-12,max_step=0.02,t_eval=xs)
    u=xs*s.y[0]
    ph=w*xs - E*E/(2*w*xs)                          # true-tail (converged) representation
    M=np.column_stack([np.sin(ph),np.cos(ph)]); (A,B),*_=np.linalg.lstsq(M,u,rcond=None)
    return np.mod(np.arctan2(B,A),np.pi)
gam=Z[:5]
print(f"X={X}, true-tail rep. Per-zero phase phi_k(beta); spread=max-min over the 5 zeros.")
print("If some beta drives spread->0, the U(2) (x=1 log knob) makes the zeros co-spectral.\n")
best=(9,None)
for beta in np.linspace(0,np.pi,13):
    phis=np.array([phase(g,beta) for g in gam])
    # spread on the circle: min over rotation of (max-min) of wrapped values
    sp=phis.max()-phis.min()
    if sp<best[0]: best=(sp,beta)
    print(f"  beta={beta:5.3f}: phi={np.round(phis,3)}  spread={sp:.4f}")
print(f"\nBEST: beta={best[1]:.3f} spread={best[0]:.4f} rad  (beta=0 Friedrichs spread was ~2.66)")
