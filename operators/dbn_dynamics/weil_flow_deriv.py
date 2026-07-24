# (b) d/ds Q_s(v) for the Weil form under the flow, fixed test v.  Q_s(v)=sum_j (2/x_j^2)|vhat'(x_j)|^2.
# Also verify exact conservation laws: d/ds sum x_j = 0 (com), d/ds sum x_j^2 = 4 N_pos (2 N_pos -1).
import numpy as np
G=[]
with open('/Users/scones/projects/riemann/data/zeros_100k.txt') as f:
    for line in f:
        s=line.strip()
        if s:
            try: G.append(float(s.split()[0]))
            except: pass
        if len(G)>=300: break
G=np.array(sorted(set(G))[:200]); X=np.concatenate([-G[::-1],G]); Npos=len(G); n=len(X)
def vel(X): D=X[:,None]-X[None,:]; np.fill_diagonal(D,np.inf); return 2*np.sum(1.0/D,axis=1)
V=vel(X)
# conservation laws
print("EXACT conservation laws (verify):")
print(f"  d/ds sum x_j = {np.sum(V):.3e}  (should be 0: center of mass conserved)")
d_m2=np.sum(2*X*V); print(f"  d/ds sum x_j^2 = {d_m2:.1f}  vs 4*Npos*(2Npos-1)={4*Npos*(2*Npos-1)}  (match => exact)")
# Weil form derivative for test v: vhat'(xi)=g(xi) (a chosen profile). Q=sum (2/x^2) g(x)^2.
# dQ/ds = sum [ (-4/x^3) g(x)^2 * v_j + (2/x^2) 2 g(x) g'(x) v_j ],  v_j=vel_j
def Qderiv(prof, dprof):
    Xp=X[X>0]; Vp=V[X>0]  # positive zeros (symmetric)
    g=prof(Xp); gp=dprof(Xp)
    dQ = np.sum((-4/Xp**3)*g**2*Vp + (2/Xp**2)*2*g*gp*Vp)
    Q  = np.sum((2/Xp**2)*g**2)
    return Q, dQ
print("\nd/ds Q_s(v) for test functions vhat'(xi)=exp(-(xi-xi0)^2/(2sig^2)) (bump probing zeros near xi0):")
sig=3.0
for xi0 in [14,18,25,35,60,100,200]:
    prof=lambda x,x0=xi0: np.exp(-(x-x0)**2/(2*sig**2))
    dprof=lambda x,x0=xi0: -(x-x0)/sig**2*np.exp(-(x-x0)**2/(2*sig**2))
    Q,dQ=Qderiv(prof,dprof)
    print(f"  xi0={xi0:3d}: Q={Q:.4e}  dQ/ds={dQ:+.4e}  {'(+)' if dQ>0 else '(-)'}")
print("\nlook: is dQ/ds definite-signed for a class of v? (=> Weil-positivity flow constraint toward Lambda)")
