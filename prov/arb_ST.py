import numpy as np
from scipy.integrate import solve_ivp
Z=np.load('operators/zeros500.npy')
lam=np.sqrt(2.0); FP=(2*np.pi*lam)**2; w=2*np.pi*lam; sp=lam; eps=1e-6
def rhs_full(alpha):
    def rhs(x,y): return [y[1],(2*x*y[1]+(FP*x*x-alpha)*y[0])/(2-x*x)]
    return rhs
def psi(x,alpha): return w*x - (16*np.pi**2 - alpha)/(2*w*x)
def Bval(gamma,X):
    alpha=-(gamma/2.0)**2; dxi0=-(FP*2-alpha)/(2*sp)
    xs=np.linspace(X-3/lam,X,400)
    R=solve_ivp(rhs_full(alpha),[sp+eps,X],[1.0,dxi0],method='DOP853',rtol=1e-12,atol=1e-13,t_eval=xs).y[0]
    u=xs*R; M=np.column_stack([np.sin(psi(xs,alpha)),np.cos(psi(xs,alpha))])
    (A,B),*_=np.linalg.lstsq(M,u,rcond=None); return B
def ladder(X):
    gs=np.linspace(12,96,680); Bs=np.array([Bval(g,X) for g in gs]); r=[]
    for i in range(len(gs)-1):
        if Bs[i]*Bs[i+1]<0:
            a,b=gs[i],gs[i+1]
            for _ in range(40):
                m=(a+b)/2
                if Bval(a,X)*Bval(m,X)<0: b=m
                else: a=m
            r.append((a+b)/2)
    return np.array(r)
def meanN(T): return (T/(2*np.pi))*np.log(T/(2*np.pi)) - T/(2*np.pi) + 7/8
def Ncount(evs,T): return np.searchsorted(evs,T,side='right')
true=Z[(Z>=12)&(Z<=96)]
Tg=np.linspace(14.5,94.5,801)  # avoid sitting exactly on a zero
S_true=np.array([Ncount(true,T) for T in Tg]) - meanN(Tg)
# mpmath argzeta S(T) cross-check on a coarse grid
import mpmath as mp; mp.mp.dps=25
Tc=np.linspace(14.5,94.5,60)
S_arg=np.array([float(mp.mpf(1)/mp.pi*mp.arg(mp.zeta(mp.mpf(1)/2+1j*mp.mpf(T)))) for T in Tc])
S_true_c=np.array([Ncount(true,T) for T in Tc]) - meanN(Tc)
print(f"[S(T) cross-check] corr(N_true-<N>, (1/pi)arg zeta) = {np.corrcoef(S_true_c,S_arg)[0,1]:.4f}  RMS={np.sqrt(np.mean((S_true_c-S_arg)**2)):.4f}")
for X in (80.0,100.0):
    evs=np.sort(ladder(X))
    F_op=np.array([Ncount(evs,T) for T in Tg]) - meanN(Tg)
    c=np.array([Ncount(evs,T)-Ncount(true,T) for T in Tg])  # integer count diff
    rho=np.corrcoef(F_op,S_true)[0,1]
    rms=np.sqrt(np.mean((F_op-S_true)**2))
    sp_mean=2*np.pi/np.log(np.mean(Tg)/(2*np.pi))
    print(f"\n=== X={X:.0f}: N_op={len(evs)} N_true={len(true)} ===")
    print(f"  rho(F_op,S) = {rho:.4f}   RMS(F_op-S)=RMS(c) = {rms:.4f}   0.5*mean-spacing = {0.5*sp_mean:.3f}")
    print(f"  max|F_op|={np.max(np.abs(F_op)):.3f}  max|S|={np.max(np.abs(S_true)):.3f}  count-diff c: mean={c.mean():+.3f} std={c.std():.3f} max|c|={np.max(np.abs(c))}")
    # also correlate against mpmath S on coarse grid
    F_c=np.array([Ncount(evs,T) for T in Tc]) - meanN(Tc)
    print(f"  rho(F_op, argzeta-S) = {np.corrcoef(F_c,S_arg)[0,1]:.4f}")
    v="A(CARRIES)" if (rho>=0.85 and rms<0.5*sp_mean) else ("B(SIZE-ONLY)" if (np.max(np.abs(F_op))<3 and rho<0.5) else "AMBIG/other")
    print(f"  VERDICT(float,this X): {v}")
