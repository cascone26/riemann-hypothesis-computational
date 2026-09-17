import RealSmoothedPNT
import DesmoothEShift
import PsiDesmoothSqueeze

set_option maxHeartbeats 4000000

open Chebyshev

/-!
# §1704 ★★★★★ — THE PNT FOR ψ WITH DE LA VALLÉE POUSSIN ERROR: `ψ(x) = x + O(x·e^{-c'√log x})`

The de-smoothing finale.  With `F = e^{-(c/2)√log x}` and window `h = x·F` (so `F² = E(x)` and
`h ≤ x/2` for large `x`), the squeeze §1609 traps `ψ(x)` between difference quotients of `Ψℝ`,
§1704a evaluates all three `Ψℝ`-errors, §1704b re-scales the shifted ones to `E(x)`, and every
term collapses to a multiple of `x²F²` — so `|ψ(x) − x| ≤ C'·x·e^{-(c/2)√log x}`.  This is the
prime number theorem for the Chebyshev function with the classical dlVP error term, EXIT0 from
Mathlib's axioms through the entire riemann9 chain.  (Report-tier honesty: formalizes KNOWN
mathematics — verified scaffolding, not progress on RH itself.)
-/

theorem pnt_psi_error :
    ∃ C c x₀ : ℝ, 0 < C ∧ 0 < c ∧ 1 ≤ x₀ ∧ ∀ x : ℝ, x₀ ≤ x →
      |ψ x - x| ≤ C * x * Real.exp (-(c * Real.sqrt (Real.log x))) := by
  obtain ⟨C, c, x₀, hC, hc, hx₀, hbound⟩ := real_smoothed_pnt
  refine ⟨1 / 2 + C * (4 + Real.exp (c * Real.sqrt (Real.log 2))), c / 2,
    max (2 * x₀) (Real.exp ((2 * Real.log 2 / c) ^ 2)),
    by positivity, by positivity,
    le_trans (by linarith) (le_max_left _ _), ?_⟩
  intro x hx
  have hx2x₀ : 2 * x₀ ≤ x := le_trans (le_max_left _ _) hx
  have hxexp : Real.exp ((2 * Real.log 2 / c) ^ 2) ≤ x := le_trans (le_max_right _ _) hx
  have hx2 : (2 : ℝ) ≤ x := by linarith
  have hx0 : (0 : ℝ) < x := by linarith
  have hx₀x : x₀ ≤ x := by linarith
  have hlog0 : (0 : ℝ) ≤ Real.log x := Real.log_nonneg (by linarith)
  have hlogge : (2 * Real.log 2 / c) ^ 2 ≤ Real.log x := (Real.le_log_iff_exp_le hx0).mpr hxexp
  have hsqge : 2 * Real.log 2 / c ≤ Real.sqrt (Real.log x) := by
    have h1 := Real.sqrt_le_sqrt hlogge
    rw [Real.sqrt_sq_eq_abs] at h1
    exact le_trans (le_abs_self _) h1
  have hF12 : Real.exp (-(c / 2 * Real.sqrt (Real.log x))) ≤ 1 / 2 := by
    have h2 := mul_le_mul_of_nonneg_left hsqge (by positivity : (0 : ℝ) ≤ c / 2)
    have h3 : c / 2 * (2 * Real.log 2 / c) = Real.log 2 := by
      field_simp
    have h1 : Real.log 2 ≤ c / 2 * Real.sqrt (Real.log x) := by
      rw [h3] at h2; exact h2
    calc Real.exp (-(c / 2 * Real.sqrt (Real.log x)))
        ≤ Real.exp (-(Real.log 2)) := Real.exp_le_exp.mpr (by linarith)
      _ = 1 / 2 := by
          rw [Real.exp_neg, Real.exp_log (by norm_num : (0 : ℝ) < 2)]
          norm_num
  set F := Real.exp (-(c / 2 * Real.sqrt (Real.log x))) with hFdef
  set Ex := Real.exp (-(c * Real.sqrt (Real.log x))) with hExdef
  set G := Real.exp (c * Real.sqrt (Real.log 2)) with hGdef
  have hF0 : (0 : ℝ) < F := by rw [hFdef]; exact Real.exp_pos _
  have hEx0 : (0 : ℝ) < Ex := by rw [hExdef]; exact Real.exp_pos _
  have hG0 : (0 : ℝ) < G := by rw [hGdef]; exact Real.exp_pos _
  have hFF : F * F = Ex := by
    rw [hFdef, hExdef, ← Real.exp_add]; congr 1; ring
  have hd0 : (0 : ℝ) < x * F := mul_pos hx0 hF0
  have hdx : x * F ≤ x / 2 := by
    have h1 := mul_le_mul_of_nonneg_left hF12 hx0.le
    linarith
  have hxmh₀ : x₀ ≤ x - x * F := by linarith
  have hxph₀ : x₀ ≤ x + x * F := by linarith
  obtain ⟨hEp, hEm⟩ := desmooth_E_shift hc hx2 hd0 hdx
  have hbx := abs_le.mp (hbound x hx₀x)
  have hbp := abs_le.mp (hbound (x + x * F) hxph₀)
  have hbm := abs_le.mp (hbound (x - x * F) hxmh₀)
  -- shifted Ψ-errors re-scaled to E(x)
  have hbp2 : psi1R (x + x * F) - (x + x * F) ^ 2 / 2 ≤ 3 * C * x ^ 2 * Ex := by
    have h1 := hbp.2
    have h2 : C * (x + x * F) ^ 2
          * Real.exp (-(c * Real.sqrt (Real.log (x + x * F))))
        ≤ C * (x + x * F) ^ 2 * Ex :=
      mul_le_mul_of_nonneg_left hEp (by positivity)
    have hsq : (x + x * F) ^ 2 ≤ 3 * x ^ 2 := by nlinarith [hdx, hx0, hd0.le]
    have h3 := mul_le_mul_of_nonneg_right hsq (mul_nonneg hC.le hEx0.le)
    nlinarith [h1, h2, h3]
  have hbm2 : psi1R (x - x * F) - (x - x * F) ^ 2 / 2 ≤ C * G * x ^ 2 * Ex := by
    have h1 := hbm.2
    have h2 : C * (x - x * F) ^ 2
          * Real.exp (-(c * Real.sqrt (Real.log (x - x * F))))
        ≤ C * (x - x * F) ^ 2 * (G * Ex) :=
      mul_le_mul_of_nonneg_left hEm (by positivity)
    have hsq : (x - x * F) ^ 2 ≤ x ^ 2 := by nlinarith [hdx, hx0, hd0.le]
    have h3 := mul_le_mul_of_nonneg_right hsq
      (mul_nonneg (mul_nonneg hC.le hG0.le) hEx0.le)
    nlinarith [h1, h2, h3]
  -- the squeeze, cleared of denominators
  have hs := psi_desmooth_squeeze (x := x) hd0
  have hs1 : psi1R x - psi1R (x - x * F) ≤ ψ x * (x * F) := (div_le_iff₀ hd0).mp hs.1
  have hs2 : ψ x * (x * F) ≤ psi1R (x + x * F) - psi1R x := (le_div_iff₀ hd0).mp hs.2
  rw [abs_le]
  constructor
  · rw [neg_le, neg_sub]
    apply le_of_mul_le_mul_right _ hd0
    nlinarith [hs1, hbm2, hbx.1, hFF, hEx0.le, hF0.le, hx0.le,
      mul_nonneg (mul_nonneg hC.le (sq_nonneg x)) hEx0.le,
      mul_nonneg (mul_nonneg (mul_nonneg hC.le hG0.le) (sq_nonneg x)) hEx0.le]
  · apply le_of_mul_le_mul_right _ hd0
    nlinarith [hs2, hbp2, hbx.1, hFF, hEx0.le, hF0.le, hx0.le,
      mul_nonneg (mul_nonneg hC.le (sq_nonneg x)) hEx0.le,
      mul_nonneg (mul_nonneg (mul_nonneg hC.le hG0.le) (sq_nonneg x)) hEx0.le]

#print axioms pnt_psi_error
