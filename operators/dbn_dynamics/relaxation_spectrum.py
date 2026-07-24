# Log-gas Laplacian L (relaxation modes of the dBN flow): L_{jk}=-2/(x_j-x_k)^2 offdiag, row-sums on diag.
# Compare ZETA zeros vs GUE vs EQUAL-spacing (same N, same mean density in a local window).
# Slowest mode (smallest nonzero eig) = long-wavelength density stiffness; distribution = relaxation spectrum.
import numpy as np
def L_of(X):
    X=np.asarray(X,float); n=len(X); D=X[:,None]-X[None,:]
    W=np.where(np.eye(n,dtype=bool),0.0,-2.0/np.where(np.eye(n,dtype=bool),1,D)**2)
    np.fill_diagonal(W,-W.sum(1)); return W
def spec(X):
    ev=np.linalg.eigvalsh(L_of(X)); ev=np.sort(ev)
    return ev[1], ev[-1], ev  # smallest-nonzero, largest, all
# ZETA: window of consecutive zeros at moderate height (near-constant density)
G=[]
with open('/Users/scones/projects/riemann/data/zeros_100k.txt') as f:
    for line in f:
        s=line.strip()
        if s:
            try: G.append(float(s.split()[0]))
            except: pass
        if len(G)>=1200: break
W=np.array(G[500:650])                 # 150 consecutive zeros around height ~ G[575]
mean_gap=np.mean(np.diff(W))
z_lo,z_hi,z_ev=spec(W)
# EQUAL spacing, same mean gap, same N, same center
c=W.mean(); n=len(W); EQ=c+(np.arange(n)-n/2)*mean_gap
e_lo,e_hi,e_ev=spec(EQ)
# GUE: eigenvalues of a random GUE matrix, scaled to same mean gap in the bulk
rng=np.random.RandomState(0); A=rng.randn(n,n)+1j*rng.randn(n,n); A=(A+A.conj().T)/2
gue=np.sort(np.linalg.eigvalsh(A)); gue=gue*(mean_gap/np.mean(np.diff(gue[n//4:3*n//4]))); gue=gue-gue.mean()+c
g_lo,g_hi,g_ev=spec(gue)
print(f"local window: {n} zeros around height {c:.0f}, mean gap {mean_gap:.4f}")
print(f"{'system':>12} {'smallest-nonzero L-eig (slow mode)':>34} {'largest L-eig (fast)':>22}")
print(f"{'ZETA':>12} {z_lo:>34.5f} {z_hi:>22.1f}")
print(f"{'equal-space':>12} {e_lo:>34.5f} {e_hi:>22.1f}")
print(f"{'GUE':>12} {g_lo:>34.5f} {g_hi:>22.1f}")
# relaxation spectral density: quantiles
print("\nL-eigenvalue quantiles (10/25/50/75/90%):")
for nm,ev in [('ZETA',z_ev),('equal',e_ev),('GUE',g_ev)]:
    q=np.percentile(ev[1:],[10,25,50,75,90])
    print(f"  {nm:>6}: {np.round(q,3)}")
print("\nzeta between equal-spacing (rigid) & GUE? The slow-mode stiffness = where criticality/Lambda lives.")
