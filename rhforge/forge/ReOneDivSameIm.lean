import Mathlib.Analysis.SpecialFunctions.Complex.Log

/-!
# §1313 — Same-height zero term is real: `Re(1/((σ+iγ)-(β+iγ))) = 1/(σ-β)`

At the DLVP test point `s = σ+iγ` and a hypothetical zero `ρ₀ = β+iγ` at the SAME height `γ`, the difference
`s-ρ₀ = σ-β` is real, so `1/(s-ρ₀) = 1/(σ-β)` is real and its real part is exactly `1/(σ-β)`.  This is the value
that turns §1312's abstract `(1/(s-ρ₀)).re` into the concrete `1/(σ-β)` appearing in `hzero` (§1311).
EXIT0 / std-axioms.
-/

open Complex

theorem re_one_div_sub_same_im (σ β γ : ℝ) :
    (1 / (((σ : ℂ) + I * (γ : ℝ)) - ((β : ℂ) + I * (γ : ℝ)))).re = 1 / (σ - β) := by
  have h : ((σ : ℂ) + I * (γ : ℝ)) - ((β : ℂ) + I * (γ : ℝ)) = ((σ - β : ℝ) : ℂ) := by
    push_cast; ring
  rw [h, ← Complex.ofReal_one, ← Complex.ofReal_div, Complex.ofReal_re]

#print axioms re_one_div_sub_same_im
