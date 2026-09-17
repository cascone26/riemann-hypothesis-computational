"""CONSTRUCTION attempt: cast the Suzuki/Weil apparatus into Srednicki class-C form H=A+σ·S and ask
the WALL-RELEVANT question: does the class-C structure change the positivity/coercivity (our λ_a wall)?
Pieces (from build_Q): arch(=pole+core+tail), prime — REAL SYMMETRIC = Srednicki S-type.
A = momentum (i d/dx) in the sine basis = imaginary-antisymmetric = Srednicki A-type. σ = Pauli (spin-½)."""
import numpy as np, sys, os
sys.path.insert(0,'operators')
from suzuki_exact_Q import build_Q
np.set_printoptions(suppress=True)

def momentum(K,a):  # <v_m|d/dx|v_n> in sine basis on [-a,a]->[0,L], real antisymmetric
    L=2*a; D=np.zeros((K,K))
    for m in range(1,K+1):
        for n in range(1,K+1):
            if m!=n and (m+n)%2==1:
                D[m-1,n-1]=(2.0*m*n/(L))*(1.0/(n-m)-1.0/(n+m))*1.0  # standard sine-basis derivative
    return (D-D.T)/2

sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]],complex); sz=np.array([[1,0],[0,-1]],complex); I2=np.eye(2,dtype=complex)

def phs_classC(H):
    """AZ class-C particle-hole: C=σy⊗I·K, check σy H* σy = -H (on the spin factor)."""
    d=H.shape[0]; n=d//2
    SY=np.kron(sy,np.eye(n,dtype=complex))
    return np.max(np.abs(SY@H.conj()@SY + H))

for a in [1.0,1.5,2.0,2.5,3.0]:
    K=max(80,int(30*a)+40)
    Q,P=build_Q(a,K,NG=400,return_pieces=True)
    arch=P["pole"]+P["core"]+P["tail"]; prime=P["prime"]
    lam_a=np.linalg.eigvalsh(Q)[0]/a
    Dm=momentum(K,a); A=1j*Dm
    s=np.linalg.norm(Q,2)
    def N(M): 
        nn=np.linalg.norm(M,2); return M/nn if nn>0 else M
    # Srednicki H = A⊗I + σx⊗arch + σz⊗prime  (arch,prime = S-type); scale-balanced
    H = np.kron(I2,A) + np.kron(sx,s*N(arch)) + np.kron(sz,s*N(prime))
    H=(H+H.conj().T)/2
    ev=np.sort(np.linalg.eigvalsh(H).real)
    gap=np.min(np.abs(ev))               # smallest |eigenvalue| = the near-0 gap
    phs=phs_classC(H)
    symm=np.max(np.abs(ev+ev[::-1]))     # spectrum ±symmetric?
    print(f"a={a:4.1f}  λ_a(Suzuki)={lam_a:+.4e}  |  class-C H: PHS-err={phs:.1e}  ±sym={symm:.1e}  near-0 gap={gap:.4e}  n_neg(H)={np.sum(ev<-1e-9)}")
print("\nKEY QUESTION: is the class-C H's near-0 gap tied to λ_a? does PHS force sign structure that λ_a lacks?")
print("If PHS-err~0, H is genuinely class C (spectrum ±sym, spin-½ present). Then the near-0 gap is the")
print("symplectic-protected quantity — compare its POSITIVITY behavior vs the raw λ_a wall.")
