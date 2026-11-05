import GOff0Symm
import DLNKOff0
import STOffFinalAssembly

open Matrix

theorem C_matrix_orthogonal_to_e_off (K : ℕ) (L : ℝ) (z : Fin K → ℂ) :
    (C_matrix (2 * (2 * K) + 1) (by omega))ᵀ.mulVec (e_off K L z) = 0 := by
  sorry

#print axioms C_matrix_orthogonal_to_e_off

theorem G_off0_isUnit_from_orthogonality (K : ℕ) (L : ℝ) (hL : L ≠ 0) (z : Fin K → ℂ) 
    (hz : ∀ k, (z k).im ≠ 0)
    (h_orth : (C_matrix (2 * (2 * K) + 1) (by omega))ᵀ.mulVec (e_off K L z) = 0) :
    IsUnit (G_off0 (2 * K) L K z).det := by
  sorry

#print axioms G_off0_isUnit_from_orthogonality

theorem G_off0_isUnit_assumed (K : ℕ) (L : ℝ) (hL : L ≠ 0) (z : Fin K → ℂ) 
    (hz : ∀ k, (z k).im ≠ 0) :
    IsUnit (G_off0 (2 * K) L K z).det := by
  sorry

#print axioms G_off0_isUnit_assumed
