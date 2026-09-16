import numpy as np

def quat_to_complex2x2(a,b,c,d):
    # q = a+bi+cj+dk  ->  standard complex 2x2 embedding
    return np.array([[a+1j*b, c+1j*d],
                      [-c+1j*d, a-1j*b]])

def embed_quat_hermitian(A,B,C,D):
    """A,B,C,D real NxN matrices representing a quaternionic matrix Q = A+Bi+Cj+Dk.
    Q is quaternionic-Hermitian iff A=A^T, and B,C,D are related to make Q_ij = conj(Q_ji):
    conj(a+bi+cj+dk)=a-bi-cj-dk, so Q_ji = conj(Q_ij) means A_ji=A_ij, B_ji=-B_ij, C_ji=-C_ij, D_ji=-D_ij
    i.e. A symmetric, B,C,D antisymmetric."""
    N = A.shape[0]
    M = np.zeros((2*N,2*N), dtype=complex)
    for i in range(N):
        for j in range(N):
            block = quat_to_complex2x2(A[i,j],B[i,j],C[i,j],D[i,j])
            M[2*i:2*i+2, 2*j:2*j+2] = block
    return M

np.random.seed(0)
N = 5
A = np.random.randn(N,N); A = (A+A.T)/2       # symmetric
B = np.random.randn(N,N); B = (B-B.T)/2       # antisymmetric
C = np.random.randn(N,N); C = (C-C.T)/2
D = np.random.randn(N,N); D = (D-D.T)/2

M = embed_quat_hermitian(A,B,C,D)
# check Hermitian
print("Hermitian check (max |M-M^dagger|):", np.max(np.abs(M - M.conj().T)))
evals = np.linalg.eigvalsh(M)
print("Eigenvalues (2N=%d of them):" % (2*N))
print(np.round(evals,10))
# check doubling: sort and compare consecutive pairs
evals_sorted = np.sort(evals)
pairs_equal = np.allclose(evals_sorted[0::2], evals_sorted[1::2], atol=1e-9)
print("Consecutive-pair doubling confirmed:", pairs_equal)
n_neg = np.sum(evals < -1e-9)
print("Number of negative eigenvalues (kappa):", n_neg, " -- even?", n_neg%2==0)

# Now try to FORCE an odd kappa by construction: take a random quaternionic Hermitian matrix,
# and check across 20 random seeds whether kappa is EVER odd.
print("\n--- stress test across 20 random quaternionic-Hermitian matrices ---")
odd_found = False
for seed in range(20):
    rng = np.random.RandomState(seed)
    N=6
    A = rng.randn(N,N); A=(A+A.T)/2
    B = rng.randn(N,N); B=(B-B.T)/2
    C = rng.randn(N,N); C=(C-C.T)/2
    D = rng.randn(N,N); D=(D-D.T)/2
    M = embed_quat_hermitian(A,B,C,D)
    ev = np.linalg.eigvalsh(M)
    k = int(np.sum(ev<-1e-9))
    if k%2!=0:
        odd_found=True
    print(f"seed={seed} kappa={k} even={k%2==0}")
print("\nAny odd kappa found across 20 random trials:", odd_found)
