import Mathlib
import GammaNormBound

/-!
# S2995 — growth of 1/Γ on the strip 0 ≤ Re z ≤ 1/2 (lane 1: hB input for three-lines)

ζ = Λ/Γℝ with Γℝ(s) = π^{−s/2}Γ(s/2); boundedness of E·G_T on the critical strip (S2984's hB) needs |1/Γ(s/2)| ≪ e^{c|t|}.
* `Gamma_le_inv` — 0 < x ≤ 1 ⇒ Γ(x) ≤ 1/x   (Γ(x+1) ≤ max(Γ1,Γ2) = 1 by convexity of Γ on [1,2]).
* `norm_sin_le_exp` — ‖sin w‖ ≤ exp|Im w|.
* `norm_inv_Gamma_le` — 0 ≤ Re z ≤ 1/2 ⇒ ‖Γ(z)⁻¹‖ ≤ (2/π)·exp(π|Im z|)
  (reflection Γ(z)Γ(1−z) = π/sin(πz), ‖Γ(1−z)‖ ≤ Γ(Re(1−z)) ≤ 1/Re(1−z) ≤ 2).
With z = s/2 (0 ≤ σ ≤ 1): |1/Γ(s/2)| ≤ (2/π)e^{π|t|/2} — the Gaussian G_T(s) = e^{(s−iT)²} kills this.
-/

open Complex

namespace InvGammaStrip

theorem Gamma_le_inv {x : ℝ} (hx0 : 0 < x) (hx1 : x ≤ 1) : Real.Gamma x ≤ 1 / x := by
  have hconv := Real.convexOn_Gamma
  have hmem : x + 1 ∈ segment ℝ (1 : ℝ) 2 := by
    rw [segment_eq_Icc (by norm_num)]; constructor <;> linarith
  have h1 : (1 : ℝ) ∈ Set.Ioi 0 := by norm_num
  have h2 : (2 : ℝ) ∈ Set.Ioi 0 := by norm_num
  have hle := hconv.le_on_segment h1 h2 hmem
  rw [Real.Gamma_one, Real.Gamma_two, max_self] at hle
  rw [Real.Gamma_add_one hx0.ne'] at hle
  rw [le_div_iff₀ hx0]; linarith

theorem norm_sin_le_exp (w : ℂ) : ‖Complex.sin w‖ ≤ Real.exp |w.im| := by
  have hsin : Complex.sin w = (Complex.exp (-w * I) - Complex.exp (w * I)) * I / 2 := rfl
  rw [hsin, norm_div, norm_mul, Complex.norm_I, mul_one, Complex.norm_two]
  have h1 : ‖Complex.exp (-w * I)‖ = Real.exp w.im := by
    rw [Complex.norm_exp]; congr 1; simp
  have h2 : ‖Complex.exp (w * I)‖ = Real.exp (-w.im) := by
    rw [Complex.norm_exp]; congr 1; simp
  have htri := norm_sub_le (Complex.exp (-w * I)) (Complex.exp (w * I))
  rw [h1, h2] at htri
  have ha : Real.exp w.im ≤ Real.exp |w.im| := Real.exp_le_exp.mpr (le_abs_self _)
  have hb : Real.exp (-w.im) ≤ Real.exp |w.im| := Real.exp_le_exp.mpr (neg_le_abs _)
  rw [div_le_iff₀ (by norm_num)]
  linarith

theorem norm_inv_Gamma_le {z : ℂ} (h0 : 0 ≤ z.re) (h1 : z.re ≤ 1 / 2) :
    ‖(Complex.Gamma z)⁻¹‖ ≤ 2 / Real.pi * Real.exp (Real.pi * |z.im|) := by
  have hpi := Real.pi_pos
  have hRHS : 0 ≤ 2 / Real.pi * Real.exp (Real.pi * |z.im|) := by positivity
  by_cases hG : Complex.Gamma z = 0
  · rw [hG, inv_zero, norm_zero]; exact hRHS
  have hre1 : 0 < (1 - z).re := by simp; linarith
  have hG1 : Complex.Gamma (1 - z) ≠ 0 := Complex.Gamma_ne_zero_of_re_pos hre1
  have hrefl := Complex.Gamma_mul_Gamma_one_sub z
  have hs : Complex.sin (↑Real.pi * z) ≠ 0 := by
    intro hs0
    rw [hs0, div_zero] at hrefl
    exact mul_ne_zero hG hG1 hrefl
  have hinv : (Complex.Gamma z)⁻¹ = Complex.Gamma (1 - z) * Complex.sin (↑Real.pi * z) / ↑Real.pi := by
    have hpiC : (↑Real.pi : ℂ) ≠ 0 := by exact_mod_cast hpi.ne'
    rw [eq_div_iff hs] at hrefl
    field_simp
    rw [show z * ↑Real.pi = ↑Real.pi * z by ring]
    linear_combination -hrefl
  rw [hinv, norm_div, norm_mul, Complex.norm_real, Real.norm_of_nonneg hpi.le]
  have hGb : ‖Complex.Gamma (1 - z)‖ ≤ 2 := by
    refine le_trans (norm_Gamma_le_Gamma_re hre1) ?_
    have hre : (1 - z).re = 1 - z.re := by simp
    rw [hre]
    refine le_trans (Gamma_le_inv (by linarith) (by linarith)) ?_
    rw [div_le_iff₀ (by linarith)]; linarith
  have hsb : ‖Complex.sin (↑Real.pi * z)‖ ≤ Real.exp (Real.pi * |z.im|) := by
    refine le_trans (norm_sin_le_exp _) (le_of_eq ?_)
    congr 1
    simp [abs_mul, abs_of_pos hpi]
  rw [div_le_iff₀ hpi]
  calc ‖Complex.Gamma (1 - z)‖ * ‖Complex.sin (↑Real.pi * z)‖
      ≤ 2 * Real.exp (Real.pi * |z.im|) := mul_le_mul hGb hsb (norm_nonneg _) (by norm_num)
    _ = 2 / Real.pi * Real.exp (Real.pi * |z.im|) * Real.pi := by field_simp

end InvGammaStrip

#print axioms InvGammaStrip.Gamma_le_inv
#print axioms InvGammaStrip.norm_sin_le_exp
#print axioms InvGammaStrip.norm_inv_Gamma_le
