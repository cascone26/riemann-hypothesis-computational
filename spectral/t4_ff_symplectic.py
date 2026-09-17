"""T4: in the PROVEN function-field case, is the structure manifestly class-C/symplectic, and does our own
apparatus (§982 conjugate-pairs) already encode it? Frobenius on H^1(C/F_q) preserves the Weil symplectic
pairing; eigenvalues come in pairs α, q/α (|α|=√q ⟺ RH_ff). Check concretely + connect to §982."""
import numpy as np

def frob_eigs_elliptic(q, a):
    # E/F_q, #E = q+1-a, Frobenius eigenvalues roots of x^2 - a x + q
    disc = a*a - 4*q
    if disc < 0:
        r = np.sqrt(-disc)/2; return np.array([a/2 + 1j*r, a/2 - 1j*r])
    r=np.sqrt(disc)/2; return np.array([a/2+r, a/2-r])

print("FUNCTION-FIELD SYMPLECTIC STRUCTURE (proven RH_ff case):")
print("Elliptic curves E/F_q: Frobenius eigenvalues α,β roots of x²−a·x+q.")
print(f"{'q':>4} {'a':>4} | {'α':>16} {'β':>16} | αβ=q? |α|=√q? β=q/ᾱ (symplectic pair)?")
for (q,a) in [(7,2),(11,-3),(13,4),(101,-6),(103,17)]:
    ev=frob_eigs_elliptic(q,a)
    al,be=ev
    prod=al*be
    absok=abs(abs(al)-np.sqrt(q))<1e-9
    # symplectic: β = q/conj(α) i.e. the pair is (α, q/α); since |α|²=q, q/α = conj(α)
    sympair=abs(be - q/np.conj(al))<1e-9 if abs(al)>0 else False
    print(f"{q:>4} {a:>4} | {al.real:7.3f}{al.imag:+.3f}i {be.real:7.3f}{be.imag:+.3f}i | {abs(prod-q)<1e-9}  {absok}  {sympair}")

print("""
ESTABLISHED (Weil, textbook): Frobenius acts on H¹(C) (dim 2g) preserving the WEIL SYMPLECTIC PAIRING
(cup product H¹×H¹→H²(1)); so Frob ∈ GSp(2g) with multiplier q. RH_ff ⟺ all eigenvalues |α|=√q ⟺
the symplectic operator's eigenvalues lie on the |·|=√q circle. Katz–Sarnak: families of curves have
SYMPLECTIC monodromy → the "symmetry class" is symplectic → AZ class C. This is the PROVEN case being class C.

4-WAY CONVERGENCE (honest — established facts + our work, a CONCEPTUAL bridge, NOT a construction/proof):
  (1) ff (proven): Frobenius symplectic via the WEIL PAIRING; eigenvalues in pairs α ↔ q/α.
  (2) Srednicki: the number-field Hamiltonian must be class C = symplectic (two-valued).
  (3) Deninger + our §936-937: the number-field obstruction is exactly a SYMPLECTIC/parity structure.
  (4) OUR §982: off-circle zeros enter as CONJUGATE PAIRS {u,ū} (= the reciprocal α↔q/α) — our apparatus
      ALREADY encodes the symplectic pairing; each pair = one negative square (κ counting).
All four name the SAME object: a symplectic pairing (ff: Weil cup product; nf: functional eq ρ↔1−ρ).

THE PRECISE CONSTRUCTION CHALLENGE (what number-field RH needs, now sharply stated):
  realize the symplectic pairing (the ff Weil-pairing analog) as a genuine self-adjoint CLASS-C Hamiltonian on
  Spec ℤ — i.e. an operator carrying the ρ↔1−ρ symplectic form intrinsically, whose spectrum is γₖ. That is what
  Deninger's missing cohomology would supply and what Srednicki's spin-½ encodes. Construction-hard, NOT walled.
""")
