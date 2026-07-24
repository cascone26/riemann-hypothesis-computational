# Screw function g_s(tau) under the dBN flow (= Weil form with cos-test). Hunt for PDE / monotone (Lyapunov) structure.
# g_s(tau) = -2 sum_j (1-cos(x_j(s) tau))/x_j(s)^2 ;  dx_j/ds = 2 sum_{k!=j} 1/(x_j-x_k).
import numpy as np, os
G=[]
with open('/Users/scones/projects/riemann/data/zeros_100k.txt') as f:
    for line in f:
        line=line.strip()
        if line:
            try: G.append(float(line.split()[0]))
            except: pass
        if len(G)>=350: break
G=np.array(sorted(set(G))[:250]); X0=np.concatenate([-G[::-1],G])
def deriv(X): D=X[:,None]-X[None,:]; np.fill_diagonal(D,np.inf); return 2*np.sum(1.0/D,axis=1)
def flow(X,dt,ns):
    X=X.copy()
    for _ in range(ns):
        k1=deriv(X);k2=deriv(X+dt/2*k1);k3=deriv(X+dt/2*k2);k4=deriv(X+dt*k3); X=X+dt/6*(k1+2*k2+2*k3+k4)
    return X
def g_s(X,tau):
    pos=X[X>0]; return -2*np.sum((1-np.cos(pos*tau))/pos**2)
taus=np.array([0.3,0.6,1.0,1.5,2.0,3.0])
ss=[-0.02,-0.01,0.0,0.01,0.02]
Xs={s:(flow(X0,np.sign(s)*0.0005,max(1,int(abs(s)/0.0005))) if s!=0 else X0) for s in ss}
print("g_s(tau) under flow:")
print("  tau  " + "  ".join(f"s={s:+.2f}" for s in ss) + "     d(g)/ds@0")
for tau in taus:
    vals=[g_s(Xs[s],tau) for s in ss]
    dgds=(vals[3]-vals[1])/0.02
    print(f"  {tau:.1f}  " + "  ".join(f"{v:+.5f}" for v in vals) + f"     {dgds:+.5f}")
# test heat-equation structure: is d(g)/ds ~ +d^2 g/dtau^2 or -g'' or c*tau*g' ? compare dg/ds to tau-derivatives of g_0.
print("\nTest PDE: dg/ds vs tau-derivatives of g_0(tau):")
taufine=np.linspace(0.05,3.5,400); g0=np.array([g_s(X0,t) for t in taufine])
dtau=taufine[1]-taufine[0]; g0pp=np.gradient(np.gradient(g0,dtau),dtau); g0p=np.gradient(g0,dtau)
for tau in [0.6,1.0,2.0]:
    i=np.argmin(abs(taufine-tau))
    Xp=flow(X0,0.0005,20); Xm=flow(X0,-0.0005,20)
    dgds=(g_s(Xp,tau)-g_s(Xm,tau))/0.02
    print(f"  tau={tau}: dg/ds={dgds:+.5f}  g''(tau)={g0pp[i]:+.4f}  tau*g'(tau)={tau*g0p[i]:+.4f}  g''+g'/tau={g0pp[i]+g0p[i]/tau:+.4f}")
print("\nlook: does dg/ds match a clean tau-operator on g? (=> screw fn satisfies a PDE under dBN flow = the new bridge)")
