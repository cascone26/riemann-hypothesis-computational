import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.Analysis.Normed.Field.Basic

/-!
# §1372 — the DLVP zero lies in the disk centred at the test point

Companion geometry to §1371: the zero `β + iγ` lies in `ball (σ+iγ) R` whenever `σ - β < R` (with `β < σ`),
since `dist(β+iγ, σ+iγ) = σ - β`.  This discharges §1355's `hρ_ball` (the kept zero is in the factorization
disk) for a disk centred at the test point `σ+iγ`.  EXIT0 / std-axioms.
-/

open Complex

theorem zero_mem_ball_test {σ γ β R : ℝ} (hβσ : β < σ) (h : σ - β < R) :
    ((β : ℂ) + I * (γ : ℝ)) ∈ Metric.ball ((σ : ℂ) + I * (γ : ℝ)) R := by
  rw [Metric.mem_ball, Complex.dist_eq]
  have hsub : ((β : ℂ) + I * (γ : ℝ)) - ((σ : ℂ) + I * (γ : ℝ)) = (((β - σ : ℝ)) : ℂ) := by
    push_cast; ring
  rw [hsub, Complex.norm_real, Real.norm_eq_abs, abs_of_neg (by linarith : (β : ℝ) - σ < 0)]
  linarith

#print axioms zero_mem_ball_test
