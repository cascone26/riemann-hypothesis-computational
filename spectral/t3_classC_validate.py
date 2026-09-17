"""T3 exp 1: validate Srednicki's class-C form H = A + σ·S produces the qualitative
signatures the Riemann Hamiltonian must have: spectrum symmetric about 0, a gap near 0
(symplectic signature), GUE-like bulk statistics. A imaginary-antisymmetric, S_i real-symmetric,
σ Pauli on an extra 2D spin space. This is the structural target before wiring in arithmetic."""
import numpy as np
np.random.seed(0)

def classC_H(N, Smats, A):
    sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]]); sz=np.array([[1,0],[0,-1]],complex)
    P=[sx,sy,sz]; I2=np.eye(2)
    H=np.kron(A,I2)+sum(np.kron(Smats[i],P[i]) for i in range(3))
    return H

def build_random_classC(N):
    B=np.random.randn(N,N); A=1j*(B-B.T)/2            # imaginary antisymmetric, hermitian
    S=[]
    for _ in range(3):
        M=np.random.randn(N,N); S.append((M+M.T)/2)    # real symmetric, hermitian
    return classC_H(N,S,A)

N=300
H=build_random_classC(N)
print("hermitian:", np.allclose(H,H.conj().T))
ev=np.sort(np.linalg.eigvalsh(H))
# class C => particle-hole symmetric spectrum about 0
sym_err=np.max(np.abs(np.sort(ev)+np.sort(-ev)[::-1]))
print(f"dim H = {H.shape[0]}, spectrum symmetric about 0? max|E+(-E)| = {sym_err:.2e}")
# gap near zero: smallest |eigenvalue| vs mean level spacing
absmin=np.min(np.abs(ev)); mean_spacing=np.mean(np.diff(ev))
print(f"gap near 0: min|E| = {absmin:.4f}, mean spacing = {mean_spacing:.4f}, ratio = {absmin/mean_spacing:.2f}")
# level spacing statistics (bulk) - GUE ratio r-statistic <r>~0.60 (GUE), 0.39 (Poisson)
d=np.diff(ev); r=np.minimum(d[:-1],d[1:])/np.maximum(d[:-1],d[1:])
print(f"level-spacing ratio <r> = {np.mean(r):.3f}  (GUE~0.60, GSE~0.68, Poisson~0.39)")
# compare to a plain GUE (class A) of same size: should have NO gap at 0
G=np.random.randn(2*N,2*N)+1j*np.random.randn(2*N,2*N); G=(G+G.conj().T)/2
evg=np.sort(np.linalg.eigvalsh(G)); dg=np.diff(evg)
rg=np.minimum(dg[:-1],dg[1:])/np.maximum(dg[:-1],dg[1:])
print(f"[control] plain GUE (class A): min|E|/spacing = {np.min(np.abs(evg))/np.mean(dg):.2f}, <r> = {np.mean(rg):.3f}")
print("=> class C should show a LARGER near-zero gap ratio than class A (the symplectic near-zero repulsion).")
