/-
TrigPolyPositivityZeta.lean (HP RH center, LEDGER S2921, 2026-11-08) — lane 1.
Summed form of S2920: for ANY cosine polynomial P(θ) = Σ_{k≤K} c_k cos(kθ) ≥ 0 and σ > 1,
  0 ≤ Σ_k c_k · Re(−ζ'/ζ)(σ + ikγ).
Generalizes §1358 (`dlvp_positivity_riemannZeta`, the 3-4-1 case). Instantiated at 5-8-4-1.
-/
import TrigPolyTermNonneg
import LogDerivZetaVonMangoldt

open Complex ArithmeticFunction LSeries
open scoped LSeries.notation

namespace TrigPolyDLVP

theorem trig_positivity_riemannZeta {σ γ : ℝ} (hσ : 1 < σ) {K : ℕ} (c : Fin (K + 1) → ℝ)
    (hP : ∀ θ : ℝ, 0 ≤ ∑ k, c k * Real.cos (k * θ)) :
    0 ≤ ∑ k : Fin (K + 1), c k *
      (-logDeriv riemannZeta ((σ : ℂ) + I * (((k : ℝ) * γ : ℝ) : ℂ))).re := by
  have hs : ∀ k : Fin (K + 1), 1 < ((σ : ℂ) + I * (((k : ℝ) * γ : ℝ) : ℂ)).re := by
    intro k
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]; simpa using hσ
  have H : ∀ k ∈ (Finset.univ : Finset (Fin (K + 1))),
      HasSum (fun n => c k * (term (↗Λ) ((σ : ℂ) + I * (((k : ℝ) * γ : ℝ) : ℂ)) n).re)
        (c k * (LSeries (↗Λ) ((σ : ℂ) + I * (((k : ℝ) * γ : ℝ) : ℂ))).re) :=
    fun k _ => (Complex.hasSum_re (LSeriesSummable_vonMangoldt (hs k)).hasSum).mul_left (c k)
  have Hsum := hasSum_sum H
  have hnn := Hsum.nonneg (fun n => trig_term_combination_nonneg (σ := σ) (γ := γ) c hP n)
  have E : ∀ k : Fin (K + 1), -logDeriv riemannZeta ((σ : ℂ) + I * (((k : ℝ) * γ : ℝ) : ℂ))
      = LSeries (↗Λ) ((σ : ℂ) + I * (((k : ℝ) * γ : ℝ) : ℂ)) := by
    intro k; rw [logDeriv_riemannZeta_eq_neg_vonMangoldt (hs k)]; ring
  simp_rw [E]
  exact hnn

/-- The 5-8-4-1 positivity for −ζ'/ζ. -/
theorem p5841_positivity_riemannZeta {σ γ : ℝ} (hσ : 1 < σ) :
    0 ≤ ∑ k : Fin 4, ![(5:ℝ), 8, 4, 1] k *
      (-logDeriv riemannZeta ((σ : ℂ) + I * (((k : ℝ) * γ : ℝ) : ℂ))).re := by
  refine trig_positivity_riemannZeta hσ ![5, 8, 4, 1] (fun θ => ?_)
  have := p5841_nonneg θ
  simp [Fin.sum_univ_four]
  linarith

end TrigPolyDLVP

#print axioms TrigPolyDLVP.trig_positivity_riemannZeta
#print axioms TrigPolyDLVP.p5841_positivity_riemannZeta
