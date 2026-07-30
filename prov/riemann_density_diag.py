import numpy as np
from scipy.integrate import solve_ivp
Z=np.load('operators/zeros500.npy'); FP=(2*np.pi)**2; w=2*np.pi
eps=1e-4; X=120.0
def phi(E,beta=0.0):
    mu=-E*E
    def rhs(x,y): return [y[1],(2*x*y[1]+(FP*x*x-mu)*y[0])/(1-x*x)]
    ICr=np.array([1.0,(mu-FP)/2.0]); ICl=np.array([np.log(eps),1.0/eps])
    IC=np.cos(beta)*ICr+np.sin(beta)*ICl
    xs=np.linspace(X-5,X,600)
    s=solve_ivp(rhs,[1+eps,X],IC,rtol=1e-10,atol=1e-12,max_step=0.02,t_eval=xs)
    u=xs*s.y[0]; ph=w*xs-E*E/(2*w*xs)
    M=np.column_stack([np.sin(ph),np.cos(ph)]); (A,B),*_=np.linalg.lstsq(M,u,rcond=None)
    return np.mod(np.arctan2(B,A),np.pi)
gam=Z[:5]
print("ORACLE density diagnostic: operator eigenvalue density vs zero density at each gamma_k")
print("  rho_op = (1/pi)|dphi/dE|   [eigenvalue crossings per unit E, extension-independent]")
print("  rho_true = (1/2pi)ln(gamma/2pi)  [Riemann-von Mangoldt local zero density]\n")
dE=0.01
for g in gam:
    # central diff of phi mod pi, unwrap the small step into (-pi/2,pi/2)
    p1=phi(g-dE); p2=phi(g+dE); d=p2-p1
    d=(d+np.pi/2)%np.pi-np.pi/2
    rho_op=abs(d/(2*dE))/np.pi
    rho_true=np.log(g/(2*np.pi))/(2*np.pi)
    print(f"  gamma={g:7.4f}  rho_op={rho_op:.4f}  rho_true={rho_true:.4f}  ratio={rho_op/rho_true:6.3f}")
print("\nIf ratio wanders far from 1 and is not constant, THIS D's Weyl law != N(T):")
print("the density-mismatch is the mechanism behind vN's 2.36-rad phase spread -> (C) intrinsic.")
