import DLVPCharFactor
import Mathlib.NumberTheory.LSeries.Dirichlet
import Mathlib.NumberTheory.DirichletCharacter.Bounds

/-!
# §1466 — termwise DLVP `3-4-1` nonnegativity for a Dirichlet character (`ζ³·L(·,χ)⁴·L(·,χ²)`)

Each L-series *term* is `a·(coeff)`, `a = Λ(n)·n^{-σ} ≥ 0`: `↗Λ` at `σ` gives `a`, `↗χ·↗Λ` at `σ+it` gives
`a·Re w`, `↗χ²·↗Λ` at `σ+2it` gives `a·Re(w²)`, `w = χ(n)·n^{-it}`, `‖w‖ ≤ 1`.  Sum `= a·(3+4Re w+Re w²) ≥ 0`
(§1465).  EXIT0 / std-axioms.
-/

open Complex ArithmeticFunction LSeries DirichletCharacter
open scoped LSeries.notation

theorem dlvp_char_term_nonneg {N : ℕ} (χ : DirichletCharacter ℂ N) {σ t : ℝ} (n : ℕ) :
    0 ≤ 3 * (term (↗Λ) ((σ : ℂ)) n).re
      + 4 * (term (↗χ * ↗Λ) ((σ : ℂ) + I * (t : ℝ)) n).re
      + (term (↗(χ ^ 2) * ↗Λ) ((σ : ℂ) + I * ((2 * t : ℝ) : ℂ)) n).re := by
  rcases eq_or_ne n 0 with hn0 | hn0
  · subst hn0; simp [LSeries.term]
  · have hn : (n : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hn0
    have hnpos : (0 : ℝ) < (n : ℝ) := by positivity
    set a : ℝ := (Λ n) * (n : ℝ) ^ (-σ) with ha_def
    have ha0 : 0 ≤ a := by
      rw [ha_def]; exact mul_nonneg vonMangoldt_nonneg (Real.rpow_pos_of_pos hnpos _).le
    set w : ℂ := (↗χ) n * (n : ℂ) ^ (-(I * (t : ℝ))) with hw_def
    have hnc_norm : ‖(n : ℂ) ^ (-(I * (t : ℝ)))‖ = 1 := by
      rw [Complex.norm_natCast_cpow_of_pos hn0.bot_lt]; simp
    have hwnorm : ‖w‖ ≤ 1 := by
      rw [hw_def, norm_mul, hnc_norm, mul_one]; exact DirichletCharacter.norm_le_one χ n
    have hreal : (n : ℂ) ^ (-(σ : ℂ)) = (((n : ℝ) ^ (-σ) : ℝ) : ℂ) := by
      rw [show ((n : ℂ)) = (((n : ℝ)) : ℂ) from by push_cast; ring,
        show (-(σ : ℂ)) = (((-σ) : ℝ) : ℂ) from by push_cast; ring,
        ← Complex.ofReal_cpow hnpos.le]
    have key : ∀ (g : ℕ → ℂ) (b : ℝ),
        (term (g * ↗Λ) ((σ : ℂ) + I * (b : ℝ)) n).re = a * (g n * (n : ℂ) ^ (-(I * (b : ℝ)))).re := by
      intro g b
      rw [LSeries.term_of_ne_zero hn0]
      have hcombine : (g * ↗Λ) n / (n : ℂ) ^ ((σ : ℂ) + I * (b : ℝ))
          = (a : ℂ) * (g n * (n : ℂ) ^ (-(I * (b : ℝ)))) := by
        show g n * (((Λ n : ℝ)) : ℂ) / (n : ℂ) ^ ((σ : ℂ) + I * (b : ℝ)) = _
        rw [Complex.cpow_add _ _ hn, ha_def]
        rw [div_eq_mul_inv, mul_inv, ← Complex.cpow_neg, ← Complex.cpow_neg, hreal]
        push_cast; ring
      rw [hcombine, Complex.re_ofReal_mul]
    have e0 : (term (↗Λ) ((σ : ℂ)) n).re = a := by
      have h1 := key 1 0
      simp only [Complex.ofReal_zero, mul_zero, add_zero, Pi.one_apply, one_mul, neg_zero,
        Complex.cpow_zero, Complex.one_re, mul_one] at h1
      exact h1
    have e1 : (term (↗χ * ↗Λ) ((σ : ℂ) + I * (t : ℝ)) n).re = a * w.re := by
      rw [key (↗χ) t, ← hw_def]
    have hchi2 : (↗(χ ^ 2)) n * (n : ℂ) ^ (-(I * ((2 * t : ℝ) : ℂ))) = w ^ 2 := by
      have hc2 : (↗(χ ^ 2)) n = ((↗χ) n) ^ 2 := by
        show (χ ^ 2) (n : ZMod N) = (χ (n : ZMod N)) ^ 2
        rw [sq, sq, MulChar.mul_apply]
      have hcpow : (n : ℂ) ^ (-(I * ((2 * t : ℝ) : ℂ))) = ((n : ℂ) ^ (-(I * (t : ℝ)))) ^ 2 := by
        rw [sq, ← Complex.cpow_add _ _ hn]; congr 1; push_cast; ring
      rw [hc2, hcpow, hw_def, mul_pow]
    have e2 : (term (↗(χ ^ 2) * ↗Λ) ((σ : ℂ) + I * ((2 * t : ℝ) : ℂ)) n).re = a * (w ^ 2).re := by
      rw [key (↗(χ ^ 2)) (2 * t), hchi2]
    rw [e0, e1, e2]
    nlinarith [ha0, dlvp_char_factor_nonneg hwnorm]

#print axioms dlvp_char_term_nonneg
