import DetDprimeCore
import DiagResolventForm

open Matrix
open scoped BigOperators

namespace CMKeyLemma

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- **CM 2025 Lemma "key" part (iii), diagonal form (capstone).**

For a diagonal operator `D = diagonal d` (CM's scaling operator, `d_j = j`), the perturbed
`D' = D − |Dξ⟩⟨η|`, and `s` off the diagonal spectrum (`d_j − s ≠ 0`), with `⟨η,ξ⟩ = 1`:

    `det(D' − s·1) = −s · det(D − s·1) · Σ_j η_j (d_j − s)⁻¹ ξ_j`.

With `d_j = j`, `η_j ≡ 1` the last factor is exactly `ξ̂(s) = Σ_j (j−s)⁻¹ ξ_j`, so **the eigenvalues of
`D''` are precisely the zeros of the finite Fourier transform of the minimal Weil eigenvector `ξ`** —
Connes–Moscovici's spectral realization, fully assembled and machine-verified, axiom-clean.
Composes `detDprime_core` (cy380) with `dotProduct_inv_diagonal_mulVec` (cy381). -/
theorem detDprime_diagonal (d η ξ : n → ℝ) (s : ℝ)
    (hg : ∀ j, d j - s ≠ 0) (hη : η ⬝ᵥ ξ = 1) :
    (diagonal d - Matrix.vecMulVec (diagonal d *ᵥ ξ) η - s • (1 : Matrix n n ℝ)).det
      = -s * (diagonal d - s • (1 : Matrix n n ℝ)).det
          * (∑ j, η j * ((d j - s)⁻¹ * ξ j)) := by
  have hdiag : diagonal d - s • (1 : Matrix n n ℝ) = diagonal (fun j => d j - s) := by
    rw [Matrix.smul_one_eq_diagonal, Matrix.diagonal_sub]
  have hA : IsUnit (diagonal d - s • (1 : Matrix n n ℝ)).det := by
    rw [hdiag, Matrix.det_diagonal, isUnit_iff_ne_zero, Finset.prod_ne_zero_iff]
    intro j _; exact hg j
  rw [detDprime_core (diagonal d) η ξ s hA hη]
  congr 1
  rw [hdiag, dotProduct_inv_diagonal_mulVec (fun j => d j - s) ξ η hg]

#print axioms detDprime_diagonal

end CMKeyLemma
