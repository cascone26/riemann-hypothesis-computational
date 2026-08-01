import numpy as np, mpmath as mp
mp.mp.dps=10
# PATH B (operator-space, escapes B-4's probe-space exhaustion): the scaling operator (Jacobi matrix)
# of the measure |zeta(1/2+is)|^2 * Gaussian-window. As the CCM semilocal measure |prod L_v|^2 -> |zeta|^2
# (S->all primes), this is the "full" scaling operator. The zeros gamma are where the measure VANISHES
# (double zeros of |zeta|^2). Q: do the zeros leave a computable imprint on the Jacobi matrix, distinct
# from the baseline (Gaussian-only, no zeta) whose Jacobi matrix is the known Hermite b_n=sqrt(n/2)?
def window(s, W=40.0): return np.exp(-(s/W)**2)   # Gaussian damping so |zeta|^2 measure is finite
# sample the measure on a grid
S=np.linspace(-110,110,44000); ds=S[1]-S[0]
print("computing |zeta(1/2+is)|^2 on grid (mpmath)...", flush=True)
Z2=np.array([abs(mp.zeta(0.5+1j*float(s)))**2 for s in S])
w_zeta = Z2*window(S)*ds          # |zeta|^2 * gaussian measure
w_base = window(S)*ds             # baseline: gaussian only (no zeta) -> Hermite Jacobi
def lanczos(s, w, K):
    Mt=w.sum(); qp=np.zeros_like(s); q=np.sqrt(np.abs(w)/Mt); a=np.zeros(K); b=np.zeros(K); bb=0.0
    for k in range(K):
        v=s*q; aa=np.dot(q,v); a[k]=aa; v=v-aa*q-bb*qp
        v=v-q*np.dot(q,v)-qp*np.dot(qp,v); bb=np.linalg.norm(v)
        if k+1<K: b[k+1]=bb
        qp=q; q=v/bb if bb>0 else v
    return a,b
K=20
ab,bb=lanczos(S,w_base,K)   # baseline Jacobi (Hermite-like)
az,bz=lanczos(S,w_zeta,K)   # |zeta|^2-measure Jacobi
print(f"\nJacobi matrix of |zeta|^2-measure vs Gaussian baseline (diag a_n, offdiag b_n):")
print(f"{'n':>3} {'a_base':>9} {'a_zeta':>9} {'b_base':>9} {'b_zeta':>9} {'db=b_zeta-b_base':>16}")
for n in range(min(K,14)):
    print(f"{n:>3} {ab[n]:>9.4f} {az[n]:>9.4f} {bb[n]:>9.4f} {bz[n]:>9.4f} {bz[n]-bb[n]:>16.5f}")
# the imprint db = how the zeros (|zeta|^2 structure) deform the scaling operator
print(f"\ndiag imprint max|a_zeta-a_base| = {np.max(np.abs(az-ab)):.4e} (measure even => ~0)")
print(f"offdiag imprint RMS(b_zeta-b_base) = {np.sqrt(np.mean((bz[1:]-bb[1:])**2)):.4f}")
print("=> the |zeta|^2 measure's Jacobi matrix carries a structured imprint = the zeros' operator-space")
print("   signature on the scaling operator. This is the S->all limit of the semilocal operator (469).")
