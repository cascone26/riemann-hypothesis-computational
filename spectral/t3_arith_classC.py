"""T3 exp 2a: fuse our arithmetic apparatus into Srednicki class-C form H=A+σ·S.
STRUCTURAL claim to check: our Weil-form pieces (arch, prime, pole) are REAL SYMMETRIC = Srednicki S-type;
our apparatus has NO imaginary-antisymmetric (A-type) piece; the natural A = momentum i d/dx (the Suzuki operator).
Then build a first fused class-C candidate and report its spectrum HONESTLY vs the real zeros (exploratory)."""
import numpy as np, sys, os
sys.path.insert(0, os.path.join(os.getcwd(),'operators'))
from suzuki_exact_Q import build_Q
G = np.array([float(x) for x in open("operators/_zeros150.txt").read().split()][:60])

a=1.5; K=90
Q, P = build_Q(a, K, NG=400, return_pieces=True)
arch = P["pole"] + P["core"] + P["tail"]   # archimedean + pole (smooth part)
prime = P["prime"]                          # prime part
# (1) structural: are they real-symmetric (S-type)?
def sym_err(M): return np.max(np.abs(M-M.T))
def imag_part(M): return np.max(np.abs(M.imag)) if np.iscomplexobj(M) else 0.0
print("STRUCTURAL CHECK (our apparatus vs Srednicki A+σ·S form):")
print(f"  arch:  real? imag={imag_part(arch):.1e}  symmetric? |M-Mᵀ|={sym_err(arch):.1e}  => {'S-type (real sym)' if sym_err(arch)<1e-9 and imag_part(arch)<1e-9 else 'other'}")
print(f"  prime: real? imag={imag_part(prime):.1e}  symmetric? |M-Mᵀ|={sym_err(prime):.1e}  => {'S-type (real sym)' if sym_err(prime)<1e-9 else 'other'}")
# our apparatus has NO A-type (imaginary antisymmetric) piece -> confirm by checking antisym part is ~0
print(f"  antisym part of full Q: |Q-Qᵀ|/2 = {np.max(np.abs((Q-Q.T)/2)):.1e}  => our apparatus is ALL S-type, NO A-type (no momentum).")
# (2) the missing A = momentum i d/dx in the Galerkin (sine) basis: d/dx sin(nπx/2a) matrix is imag-antisym
n = np.arange(1,K+1)
# <sin(mπx/L)| d/dx |sin(nπx/L)> on [0,L], L=2a: = (2/L)∫ sin(mπx/L)(nπ/L)cos(nπx/L)dx = coupling, antisymmetric
L=2*a
D = np.zeros((K,K))
for i in range(K):
    for j in range(K):
        if (i+1+j+1)%2==1 and i!=j:  # odd m+n -> nonzero
            D[i,j] = (2.0/L)*((j+1)*np.pi/L)*( (1-np.cos((i+1+j+1)*np.pi))/((i+1+j+1)*np.pi/L) + (1-np.cos((i+1-(j+1))*np.pi))/(((i+1)-(j+1))*np.pi/L) ) if i!=j else 0
D = (D - D.T)/2  # antisymmetrize -> real antisymmetric; A = i*D is imaginary-antisym hermitian
A = 1j*D
print(f"\n  A = i·(momentum) : imaginary-antisymmetric hermitian? |A-A†|={np.max(np.abs(A-A.conj().T)):.1e} (0=OK)")

# (3) first fused class-C candidate: H = A⊗I2 + σz⊗arch + σx⊗prime  (arch,prime as S-components)
sx=np.array([[0,1],[1,0]],complex); sz=np.array([[1,0],[0,-1]],complex); I2=np.eye(2)
# normalize scales so pieces are comparable
def nrm(M): return M/np.linalg.norm(M,2)
scale = np.linalg.norm(A,2)
H = np.kron(A,I2) + np.kron(scale*nrm(arch),sz) + np.kron(scale*nrm(prime),sx)
print(f"  H hermitian? {np.allclose(H,H.conj().T)}   dim={H.shape[0]}")
ev=np.sort(np.linalg.eigvalsh(H))
print(f"  spectrum symmetric about 0 (class-C PHS)? max|E+(-E)|={np.max(np.abs(ev+ev[::-1])):.2e}")
pos=ev[ev>1e-9]
print(f"  # positive eigenvalues: {len(pos)}; smallest few: {np.round(pos[:6],3)}")
print(f"  real γ (first 6): {np.round(G[:6],3)}")
print("\nHONEST READ: this is EXPLORATORY. No principled reason a naive scale-matched combination gives γₖ")
print("(that would BE RH). Expected: spectrum is class-C-symmetric but does NOT match γₖ. Value = the STRUCTURAL")
print("bridge: our Weil apparatus IS the S-type; the missing ingredients are A=momentum (Suzuki i d/dx) + spin-½.")
print("The principled T3 build = fuse the SUZUKI operator (whose scalar spectrum≈γₖ) with the spin, not a scale-fit.")
