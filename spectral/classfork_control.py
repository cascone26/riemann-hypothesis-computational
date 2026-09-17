"""DECISIVE class-A vs class-C control: is ζ's <r>=0.615 excess over GUE(0.5996) a real
symplectic signal, or finite-N GUE noise? Draw actual GUE (class A) AND GSE (class C) matrices
at matched N, compute the SAME unfolding-free <r> statistic, and locate ζ between them with error bars."""
import numpy as np
rng = np.random.default_rng(20260822)

z = np.load('operators/zeros10k.npy')
d = np.diff(z); r = np.minimum(d[1:],d[:-1])/np.maximum(d[1:],d[:-1])
zeta_r = r.mean()
print(f"ζ zeros (N=10000):  <r> = {zeta_r:.4f}")
# windowed trend (low->high γ)
for lo,hi in [(0,1000),(1000,3000),(3000,6000),(6000,9998)]:
    print(f"    γ-window [{lo}:{hi}]  <r> = {r[lo:hi].mean():.4f}")

def r_of_eigs(ev):
    ev = np.sort(ev); s = np.diff(ev)
    rr = np.minimum(s[1:],s[:-1])/np.maximum(s[1:],s[:-1])
    return rr.mean()

def gue_matrix(n):
    A = (rng.standard_normal((n,n)) + 1j*rng.standard_normal((n,n)))/np.sqrt(2)
    return (A+A.conj().T)/2
def gse_matrix(n):  # class C / symplectic: quaternion self-dual, use 2n x 2n symplectic
    # GSE via quaternion real structure: H = [[A, B],[-B*, A*]] hermitian, A herm, B antisym
    A = (rng.standard_normal((n,n))+1j*rng.standard_normal((n,n)))/np.sqrt(2); A=(A+A.conj().T)/2
    B = (rng.standard_normal((n,n))+1j*rng.standard_normal((n,n)))/np.sqrt(2); B=(B-B.T)/2
    top=np.hstack([A,B]); bot=np.hstack([-B.conj(),A.conj()])
    return np.vstack([top,bot])

# match the *number of spacings* to a few hundred per matrix, many trials for error bars
N=400; T=60
gue_rs=[]; gse_rs=[]
for _ in range(T):
    ev=np.linalg.eigvalsh(gue_matrix(N)); mid=ev[N//4:3*N//4]; gue_rs.append(r_of_eigs(mid))
    ev=np.linalg.eigvalsh(gse_matrix(N)); mid=ev[N//2:3*N//2]; gse_rs.append(r_of_eigs(mid))
gue_rs=np.array(gue_rs); gse_rs=np.array(gse_rs)
print(f"\nFinite-N controls (N={N}, {T} trials, bulk eigenvalues only):")
print(f"  GUE  (class A): <r> = {gue_rs.mean():.4f} ± {gue_rs.std():.4f}")
print(f"  GSE  (class C): <r> = {gse_rs.mean():.4f} ± {gse_rs.std():.4f}")
# z-score of zeta vs each (using ensemble std as the finite-N spread, plus zeta's own SE)
zeta_se = r.std()/np.sqrt(len(r))
print(f"\n  ζ <r> = {zeta_r:.4f}  (SE={zeta_se:.4f}, from {len(r)} spacings)")
print(f"  ζ vs GUE:  Δ = {zeta_r-gue_rs.mean():+.4f}   z ≈ {(zeta_r-gue_rs.mean())/np.sqrt(gue_rs.std()**2+zeta_se**2):.1f}σ")
print(f"  ζ vs GSE:  Δ = {zeta_r-gse_rs.mean():+.4f}   z ≈ {(zeta_r-gse_rs.mean())/np.sqrt(gse_rs.std()**2+zeta_se**2):.1f}σ")
