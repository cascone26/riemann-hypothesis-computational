import numpy as np
from scipy.integrate import solve_ivp
Z=np.load('operators/zeros500.npy'); FP=(2*np.pi)**2; w=2*np.pi; eps=1e-4
def phi(E,X):
    mu=-E*E
    def rhs(x,y): return [y[1],(2*x*y[1]+(FP*x*x-mu)*y[0])/(1-x*x)]
    IC=np.array([1.0,(mu-FP)/2.0])
    xs=np.linspace(X-5,X,600)
    s=solve_ivp(rhs,[1+eps,X],IC,rtol=1e-10,atol=1e-12,max_step=0.02,t_eval=xs)
    u=xs*s.y[0]; ph=w*xs-E*E/(2*w*xs)
    M=np.column_stack([np.sin(ph),np.cos(ph)]); (A,B),*_=np.linalg.lstsq(M,u,rcond=None)
    return np.mod(np.arctan2(B,A),np.pi)
gam=Z[:5]; dE=0.01
print("X-stability: is the abs density box-scaled (~X) while the E-SHAPE (ratio-drift) is stable?\n")
for X in [60.0,120.0,180.0]:
    rr=[]
    for g in gam:
        d=phi(g+dE,X)-phi(g-dE,X); d=(d+np.pi/2)%np.pi-np.pi/2
        rho_op=abs(d/(2*dE))/np.pi; rho_true=np.log(g/(2*np.pi))/(2*np.pi)
        rr.append(rho_op/rho_true)
    rr=np.array(rr)
    print(f"  X={X:5.0f}: ratio(gamma_k)={np.round(rr,3)}  drift(r1/r5)={rr[0]/rr[-1]:.3f}")
print("\nIf abs ratio grows with X but drift(r1/r5) stays ~const -> shape is intrinsic (E-form wrong),")
print("not a scale a single build constant can rescale away.")
