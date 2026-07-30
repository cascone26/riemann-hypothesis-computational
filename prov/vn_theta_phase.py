import numpy as np
from scipy.integrate import solve_ivp
Z=np.load('operators/zeros500.npy'); FP=(2*np.pi)**2; w=2*np.pi
def AB(E,X,corrected):
    mu=-E*E
    def rhs(x,y): return [y[1],(2*x*y[1]+(FP*x*x-mu)*y[0])/(1-x*x)]
    xs=np.linspace(X-5,X,600)
    s=solve_ivp(rhs,[1+1e-4,X],[1.0,(mu-FP)/2.0],rtol=1e-10,atol=1e-12,max_step=0.02,t_eval=xs)
    u=xs*s.y[0]
    ph=w*xs - (E*E/(2*w*xs) if corrected else 0.0)
    M=np.column_stack([np.sin(ph),np.cos(ph)])
    (A,B),*_=np.linalg.lstsq(M,u,rcond=None); return A,B
# phase theta = atan2(B,A) is the extension parameter the eigenfunction "wants" at energy E.
# A genuine self-adjoint extension => theta CONSTANT across the true zeros.
gam=Z[:5]  # 14.135,21.022,25.011,30.425,32.935
for X in [45,90,180]:
    print(f"\n--- X={X} ---   phase phi(E)=atan2(B,A) at the five TRUE zeros (mod pi):")
    for corr in [False,True]:
        phis=[]
        for g in gam:
            A,B=AB(g,X,corr); phis.append(np.mod(np.arctan2(B,A),np.pi))
        phis=np.array(phis)
        spread=phis.max()-phis.min()
        tag='TRUE-TAIL' if corr else 'BARE     '
        print(f"  {tag}: phi={np.round(phis,4)}  SPREAD={spread:.4f} rad  (=0 => single self-adjoint theta)")
