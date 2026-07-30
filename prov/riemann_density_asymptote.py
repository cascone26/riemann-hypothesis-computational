import numpy as np
from scipy.integrate import solve_ivp
Z=np.load('operators/zeros500.npy'); FP=(2*np.pi)**2; w=2*np.pi; eps=1e-4
def phi(E,X):
    mu=-E*E
    def rhs(x,y): return [y[1],(2*x*y[1]+(FP*x*x-mu)*y[0])/(1-x*x)]
    IC=np.array([1.0,(mu-FP)/2.0]); xs=np.linspace(X-5,X,900)
    s=solve_ivp(rhs,[1+eps,X],IC,rtol=1e-11,atol=1e-13,max_step=0.01,t_eval=xs)
    u=xs*s.y[0]; ph=w*xs-E*E/(2*w*xs)
    M=np.column_stack([np.sin(ph),np.cos(ph)]); (A,B),*_=np.linalg.lstsq(M,u,rcond=None)
    return np.mod(np.arctan2(B,A),np.pi)
idx=[0,4,9,19,49,99]; dE=0.01; X=250.0
print(f"ASYMPTOTIC density-ratio: does dphi/dE : N'(T) PEEL (->intrinsic law defect) or FLATTEN")
print(f"(->correct Weyl law, floor=low-E fluctuation, REOPENS curability)?  X={X}\n")
print(" idx   gamma     rho_op    N'(T)    ratio")
prev=None; rr=[]
for i in idx:
    g=Z[i]
    d=phi(g+dE,X)-phi(g-dE,X); d=(d+np.pi/2)%np.pi-np.pi/2
    rho=abs(d/(2*dE))/np.pi; nt=np.log(g/(2*np.pi))/(2*np.pi); r=rho/nt; rr.append(r)
    print(f" #{i:<3d} {g:8.3f}  {rho:.4f}  {nt:.4f}  {r:6.3f}")
rr=np.array(rr)
print(f"\ntrend ratio(#0->#99): {rr[0]:.3f} -> {rr[-1]:.3f}")
print("MONOTONE-FALLING toward 1  => density law approaching N(T): NOT intrinsic, floor is low-E.")
print("TURNS UP / diverges         => power-law vs log: intrinsic density defect, measure-deform forced.")
