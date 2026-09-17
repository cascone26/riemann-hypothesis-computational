import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.Analysis.Meromorphic.Divisor

/-!
# §1395 — drop-zeros upper bound: `Re(-F) ≤ -Re g` for `F = ∑ m/(z-ρ) + g`

If the partial-fractions sum has nonnegative real part termwise (`m ρ ≥ 0`, `Re(1/(z-ρ)) ≥ 0` — true at a point `z`
with `Re z > Re ρ` for every zero `ρ`), then dropping the entire sum gives `Re(-F) ≤ -Re g`.  This is the DLVP
"drop the zeros" step used for the growth bound `hgrowth`: at `z = σ+2iγ` (`σ > 1 ≥ Re ρ`) every `Re(1/(z-ρ)) > 0`,
so `Re(-logDeriv ζ(z)) ≤ -Re(logDeriv g_full z) ≤ B` (the BC bound), giving `T₂ ≤ C = B ≍ log|t|`.  EXIT0 / std.
-/

open Complex

theorem re_neg_decomp_le {S : Finset ℂ} {m : ℂ → ℤ} {z g F : ℂ}
    (hF : F = (∑ ρ ∈ S, (m ρ : ℂ) / (z - ρ)) + g)
    (hm : ∀ ρ ∈ S, 0 ≤ m ρ) (hre : ∀ ρ ∈ S, 0 ≤ (1 / (z - ρ)).re) :
    (-F).re ≤ -g.re := by
  have hsum : 0 ≤ (∑ ρ ∈ S, (m ρ : ℂ) / (z - ρ)).re := by
    rw [Complex.re_sum]
    apply Finset.sum_nonneg
    intro ρ hρ
    have hterm : ((m ρ : ℂ) / (z - ρ)).re = (m ρ : ℝ) * (1 / (z - ρ)).re := by
      rw [div_eq_mul_inv, ← one_div,
        show ((m ρ : ℂ)) = (((m ρ : ℝ)) : ℂ) from by push_cast; ring, Complex.re_ofReal_mul]
    rw [hterm]
    exact mul_nonneg (by exact_mod_cast hm ρ hρ) (hre ρ hρ)
  rw [hF, Complex.neg_re, Complex.add_re]
  linarith

#print axioms re_neg_decomp_le
