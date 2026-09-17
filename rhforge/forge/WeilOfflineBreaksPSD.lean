import Mathlib
import NotPSDIsotropicNeg

open Matrix

/-- **Backward Weil criterion (one off-line zero destroys positivity), assembled.**

`G` is the Gram sum over the on-line zeros (PSD by `WeilFormPSD`); one off-line zero contributes the
rank-2 term `vecMulVec a a - vecMulVec b b` (`a = Re w`, `b = Im w`). If a direction `x` lies in the
Gram null cone (`xᵀG x = 0`), annihilates `a` (`a⬝x = 0`) but not `b` (`b⬝x ≠ 0`) — the generic
off-line configuration — then the full Weil form `G + (aaᵀ - bbᵀ)` is NOT positive semidefinite,
since `xᵀ(G + aaᵀ - bbᵀ)x = 0 + (a⬝x)² - (b⬝x)² = -(b⬝x)² < 0`.

This is the Lean assembly of the numerical degree-2 ghost-zero detector (cycle 253): the on-line part
is harmless on `x`, and the single off-line zero injects the fatal negative direction. -/
theorem weil_offline_breaks_psd {n : Type*} [Fintype n] [DecidableEq n]
    (G : Matrix n n ℝ) (a b x : n → ℝ)
    (hGiso : x ⬝ᵥ (G *ᵥ x) = 0) (h1 : a ⬝ᵥ x = 0) (h2 : b ⬝ᵥ x ≠ 0) :
    ¬ (G + (vecMulVec a a - vecMulVec b b)).PosSemidef := by
  apply not_posSemidef_of_isotropic_neg G _ x hGiso
  have key : ∀ v : n → ℝ, x ⬝ᵥ (vecMulVec v v *ᵥ x) = (v ⬝ᵥ x) * (v ⬝ᵥ x) := by
    intro v
    simp only [dotProduct, mulVec, vecMulVec_apply]
    rw [Finset.sum_mul_sum]
    simp only [Finset.mul_sum]
    apply Finset.sum_congr rfl; intro i _
    apply Finset.sum_congr rfl; intro j _
    ring
  rw [sub_mulVec, dotProduct_sub, key a, key b, h1]
  nlinarith [mul_self_pos.mpr h2]
