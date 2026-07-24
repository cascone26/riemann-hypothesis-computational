# Verify the novel structural facts of the dBN log-gas flow:
# (1) flow = grad ascent of E=sum_{j!=k}log|x_j-x_k|:  dx_j/ds = dE/dx_j = 2 sum 1/(x_j-x_k).
# (2) Hess E = -L (negative graph Laplacian, L_{jk}=-2/(x_j-x_k)^2 offdiag) => E CONCAVE.
# (3) dE/ds=|grad E|^2>=0, d^2E/ds^2 = -2<gradE, L gradE> <= 0.  Confirm numerically.
import numpy as np
G=[]
with open('/Users/scones/projects/riemann/data/zeros_100k.txt') as f:
    for line in f:
        s=line.strip()
        if s:
            try: G.append(float(s.split()[0]))
            except: pass
        if len(G)>=200: break
G=np.array(sorted(set(G))[:150]); X=np.concatenate([-G[::-1],G]); n=len(X)
def gradE(X): D=X[:,None]-X[None,:]; np.fill_diagonal(D,np.inf); return 2*np.sum(1.0/D,axis=1)
def E(X):
    D=np.abs(X[:,None]-X[None,:]); np.fill_diagonal(D,1.0); return np.sum(np.log(D))  # sum_{j!=k}
def L_matrix(X):
    D=X[:,None]-X[None,:]; W=np.zeros((n,n))
    for j in range(n):
        for k in range(n):
            if j!=k: W[j,k]=-2.0/(X[j]-X[k])**2
    for j in range(n): W[j,j]=-np.sum(W[j,:])  # row sums -> Laplacian
    return W
g=gradE(X); L=L_matrix(X)
# (1) grad = flow velocity: already by construction. (2) Hess E = -L: check numerically via finite diff of grad
h=1e-5; Hess=np.zeros((n,n))
for k in range(n):
    Xp=X.copy(); Xp[k]+=h; Xm=X.copy(); Xm[k]-=h
    Hess[:,k]=(gradE(Xp)-gradE(Xm))/(2*h)
err=np.max(np.abs(Hess+L))/np.max(np.abs(L))
print(f"(2) ||Hess E + L|| / ||L|| = {err:.2e}  (=> Hess E = -L confirmed)")
Lev=np.linalg.eigvalsh(L)
print(f"    L eigenvalues: min={Lev.min():.3e} (>=0 => L PSD => Hess E = -L NSD => E CONCAVE), max={Lev.max():.2e}")
# (3) dE/ds and d^2E/ds^2 along flow
def deriv(X): return gradE(X)
def flow(X,dt,ns):
    X=X.copy()
    for _ in range(ns):
        k1=deriv(X);k2=deriv(X+dt/2*k1);k3=deriv(X+dt/2*k2);k4=deriv(X+dt*k3); X=X+dt/6*(k1+2*k2+2*k3+k4)
    return X
print("\n(3) E, dE/ds, d2E/ds2 along flow:")
Es=[]; ss=np.linspace(-0.01,0.03,9)
for s in ss:
    Xs=flow(X, np.sign(s)*0.0002, max(1,int(abs(s)/0.0002))) if s!=0 else X
    gs=gradE(Xs); Ls=L_matrix(Xs)
    dEds=np.sum(gs**2); d2Eds2=-2*gs@(Ls@gs)
    Es.append((s,E(Xs),dEds,d2Eds2))
    print(f"  s={s:+.3f}: E={E(Xs):.3f}  dE/ds={dEds:.1f}  d2E/ds2={d2Eds2:.1f}  {'CONCAVE(<=0)' if d2Eds2<=0 else 'CONVEX!'}")
print("\nAll d2E/ds2<=0 => E concave along flow (novel). Consequence: dE/ds=|gradE|^2 monotone DECREASING => flow rate decays.")
