import Mathlib

open Matrix

/-- **Backward Weil criterion, quantitative form (drops the exact-isotropy hypothesis).**

The realistic condition the ghost-zero detector actually uses: the off-line zero destroys positivity
whenever its negative contribution `(b⬝x)²` locally dominates the on-line part plus the off-line
positive part, `xᵀG x + (a⬝x)²`. If for some direction `x`,
`xᵀG x < (b⬝x)² - (a⬝x)²`, then the full Weil form `G + (vecMulVec a a - vecMulVec b b)` is NOT
positive semidefinite: `xᵀ(G + aaᵀ - bbᵀ)x = xᵀG x + (a⬝x)² - (b⬝x)² < 0`.

Generalizes `weil_offline_breaks_psd`: the exact-isotropy case `xᵀG x = 0` is recovered by taking
`b⬝x ≠ 0, a⬝x = 0`. No idealized null-space assumption — only that the off-line negativity outweighs
the local positivity, exactly what happens numerically once a zero moves far enough off the line. -/
theorem weil_offline_breaks_psd_quant {n : Type*} [Fintype n] [DecidableEq n]
    (G : Matrix n n ℝ) (a b x : n → ℝ)
    (hquant : x ⬝ᵥ (G *ᵥ x) < (b ⬝ᵥ x) * (b ⬝ᵥ x) - (a ⬝ᵥ x) * (a ⬝ᵥ x)) :
    ¬ (G + (vecMulVec a a - vecMulVec b b)).PosSemidef := by
  intro hpsd
  have hx := hpsd.dotProduct_mulVec_nonneg x
  rw [add_mulVec, dotProduct_add] at hx
  simp only [star_trivial] at hx
  have key : ∀ v : n → ℝ, x ⬝ᵥ (vecMulVec v v *ᵥ x) = (v ⬝ᵥ x) * (v ⬝ᵥ x) := by
    intro v
    simp only [dotProduct, mulVec, vecMulVec_apply]
    rw [Finset.sum_mul_sum]
    simp only [Finset.mul_sum]
    apply Finset.sum_congr rfl; intro i _
    apply Finset.sum_congr rfl; intro j _
    ring
  rw [sub_mulVec, dotProduct_sub, key a, key b] at hx
  linarith [hquant]
