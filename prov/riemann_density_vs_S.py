import numpy as np
from scipy.integrate import solve_ivp
Z=np.load('operators/zeros500.npy'); w0=2*np.pi; eps=1e-4; X=120.0
def phi(E,s):
    FP=(s*w0)**2; w=s*w0; mu=-E*E
    def rhs(x,y): return [y[1],(2*x*y[1]+(FP*x*x-mu)*y[0])/(1-x*x)]
    IC=np.array([1.0,(mu-FP)/2.0]); xs=np.linspace(X-5,X,600)
    sol=solve_ivp(rhs,[1+eps,X],IC,rtol=1e-10,atol=1e-12,max_step=0.02,t_eval=xs)
    u=xs*sol.y[0]; ph=w*xs-E*E/(2*w*xs)
    M=np.column_stack([np.sin(ph),np.cos(ph)]); (A,B),*_=np.linalg.lstsq(M,u,rcond=None)
    return np.mod(np.arctan2(B,A),np.pi)
gam=Z[:5]; dE=0.01
print("Does the S-minimizing frequency also flatten the density-ratio dphi/dE : N'(T)?")
print("If ratio-drift persists at the S-optimum, S-min != density-match -> shape scan must target density.\n")
for s in [1.0,1.414,1.540]:
    rr=[]
    for g in gam:
        d=phi(g+dE,s)-phi(g-dE,s); d=(d+np.pi/2)%np.pi-np.pi/2
        rho_op=abs(d/(2*dE))/np.pi; rho_true=np.log(g/(2*np.pi))/(2*np.pi)
        rr.append(rho_op/rho_true)
    rr=np.array(rr)
    print(f"  s={s:5.3f}: ratio={np.round(rr,3)}  drift(r1/r5)={rr[0]/rr[-1]:.3f}  mean={rr.mean():.3f}")
print("\nInterpretation: a scalar s rescales the ratio LEVEL (mean); it cannot remove the DRIFT")
print("(r1/r5). Drift = the E-shape mismatch = what only the x^2-weight exponent can touch.")
