/-
TrigPolyTermNonneg.lean (HP RH center, LEDGER S2920, 2026-11-08) — lane 1, beyond the 3-4-1 polynomial.
Generalizes §1357 (`dlvp_term_combination_nonneg`, the 3-4-1 case) to an ARBITRARY nonnegative cosine
polynomial P(θ) = Σ_k c_k cos(kθ): the P-combination of the von Mangoldt L-series terms at σ + ikγ is ≥ 0.
This is the input format for improved zero-free-region constants (Stečkin / Kadiri / Mossinghoff–Trudgian
use higher-degree P with larger c₁/c₀). Also proves the degree-3 P = 5 + 8cosθ + 4cos2θ + cos3θ ≥ 0
(= (1+cosθ)(1+2cosθ)², ratio c₁/c₀ = 8/5 > 4/3 of 3-4-1).
-/
import Mathlib.NumberTheory.LSeries.Dirichlet

open Complex ArithmeticFunction LSeries
open scoped LSeries.notation

namespace TrigPolyDLVP

lemma re_natCast_cpow_neg_I_mul {n : ℕ} (hn0 : n ≠ 0) (b : ℝ) :
    ((n : ℂ) ^ (-(I * (b : ℂ)))).re = Real.cos (b * Real.log n) := by
  have hn : (n : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hn0
  have hl : Complex.log (n : ℂ) = ((Real.log n : ℝ) : ℂ) := by
    rw [← Complex.ofReal_natCast, ← Complex.ofReal_log (Nat.cast_nonneg n)]
  have him : (Complex.log (n : ℂ)).im = 0 := by rw [hl]; exact Complex.ofReal_im _
  have hre : (Complex.log (n : ℂ)).re = Real.log n := by rw [hl]; exact Complex.ofReal_re _
  rw [Complex.cpow_def_of_ne_zero hn, Complex.exp_re]
  simp [Complex.mul_re, Complex.mul_im, Real.cos_neg]
  rw [him, hre]; simp [mul_comm]

lemma term_re {n : ℕ} (hn0 : n ≠ 0) (σ b : ℝ) :
    (term (↗Λ) ((σ : ℂ) + I * (b : ℂ)) n).re
      = (Λ n * (n : ℝ) ^ (-σ)) * Real.cos (b * Real.log n) := by
  have hn : (n : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hn0
  have hnpos : (0 : ℝ) < (n : ℝ) := by positivity
  rw [LSeries.term_of_ne_zero hn0, ← re_natCast_cpow_neg_I_mul hn0 b]
  have hsplit : (n : ℂ) ^ ((σ : ℂ) + I * (b : ℂ))
      = (n : ℂ) ^ (σ : ℂ) * (n : ℂ) ^ (I * (b : ℂ)) := Complex.cpow_add _ _ hn
  have hreal : (n : ℂ) ^ (-(σ : ℂ)) = (((n : ℝ) ^ (-σ) : ℝ) : ℂ) := by
    rw [show ((n : ℂ)) = (((n : ℝ)) : ℂ) from by push_cast; ring,
      show (-(σ : ℂ)) = (((-σ) : ℝ) : ℂ) from by push_cast; ring,
      ← Complex.ofReal_cpow hnpos.le]
  have hcombine : ((Λ n : ℝ) : ℂ) / ((n : ℂ) ^ (σ : ℂ) * (n : ℂ) ^ (I * (b : ℂ)))
      = (((Λ n * (n : ℝ) ^ (-σ)) : ℝ) : ℂ) * (n : ℂ) ^ (-(I * (b : ℂ))) := by
    rw [div_eq_mul_inv, mul_inv, ← Complex.cpow_neg, ← Complex.cpow_neg, hreal]
    push_cast; ring
  rw [hsplit, hcombine, Complex.re_ofReal_mul]

/-- **General nonnegative-trig-polynomial termwise positivity.** -/
theorem trig_term_combination_nonneg {σ γ : ℝ} {K : ℕ} (c : Fin (K + 1) → ℝ)
    (hP : ∀ θ : ℝ, 0 ≤ ∑ k, c k * Real.cos (k * θ)) (n : ℕ) :
    0 ≤ ∑ k : Fin (K + 1), c k * (term (↗Λ) ((σ : ℂ) + I * (((k : ℝ) * γ : ℝ) : ℂ)) n).re := by
  rcases eq_or_ne n 0 with hn0 | hn0
  · subst hn0; simp [LSeries.term]
  · have ha : 0 ≤ Λ n * (n : ℝ) ^ (-σ) :=
      mul_nonneg vonMangoldt_nonneg (Real.rpow_nonneg (Nat.cast_nonneg _) _)
    simp_rw [term_re hn0]
    have : ∑ k : Fin (K + 1), c k * ((Λ n * (n : ℝ) ^ (-σ)) * Real.cos ((k : ℝ) * γ * Real.log n))
        = (Λ n * (n : ℝ) ^ (-σ)) * ∑ k : Fin (K + 1), c k * Real.cos (k * (γ * Real.log n)) := by
      rw [Finset.mul_sum]; refine Finset.sum_congr rfl fun k _ => ?_; rw [mul_assoc]; ring_nf
    rw [this]
    exact mul_nonneg ha (hP _)

/-- The degree-3 polynomial 5 + 8cosθ + 4cos2θ + cos3θ = (1+cosθ)(1+2cosθ)² ≥ 0. -/
theorem p5841_nonneg (θ : ℝ) :
    0 ≤ 5 + 8 * Real.cos θ + 4 * Real.cos (2 * θ) + Real.cos (3 * θ) := by
  rw [Real.cos_two_mul, Real.cos_three_mul]
  have h1 : 0 ≤ 1 + Real.cos θ := by linarith [Real.neg_one_le_cos θ]
  have : 5 + 8 * Real.cos θ + 4 * (2 * Real.cos θ ^ 2 - 1) + (4 * Real.cos θ ^ 3 - 3 * Real.cos θ)
      = (1 + Real.cos θ) * (1 + 2 * Real.cos θ) ^ 2 := by ring
  rw [this]; positivity

/-- The 5-8-4-1 von Mangoldt combination is termwise nonnegative. -/
theorem p5841_term_nonneg {σ γ : ℝ} (n : ℕ) :
    0 ≤ ∑ k : Fin 4, ![5, 8, 4, 1] k *
      (term (↗Λ) ((σ : ℂ) + I * (((k : ℝ) * γ : ℝ) : ℂ)) n).re := by
  refine trig_term_combination_nonneg ![5, 8, 4, 1] (fun θ => ?_) n
  have := p5841_nonneg θ
  simp [Fin.sum_univ_four]
  linarith

end TrigPolyDLVP

#print axioms TrigPolyDLVP.trig_term_combination_nonneg
#print axioms TrigPolyDLVP.p5841_nonneg
#print axioms TrigPolyDLVP.p5841_term_nonneg
