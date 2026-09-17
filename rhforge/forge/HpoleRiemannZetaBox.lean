import ReInvLeInvRe
import LogDerivTrivChar1OneEq
import PoleConstBoundBox

set_option maxHeartbeats 1000000

/-!
# §1572 — COMPLEX `hpole`: `Re(-logDeriv ζ s) ≤ 1/(Re s - 1) + A₀` on the strip `1 < Re ≤ 2`, `|Im| ≤ 1`

Assembles the small-|t| pole bound: §1570 (`-logDeriv ζ s = 1/(s-1) - logDeriv(LFunctionTrivChar₁ 1) s`) +
§1569 (`Re(1/(s-1)) ≤ 1/(Re s-1)`) + §1571 (uniform `A₀` on the rectangle).  Unlike §1405 (real axis only), this
holds at EVERY point of the strip — in particular at the DLVP test points `σK+iγ`, `σK+2iγ` with `|γ|` SMALL,
where the T₂ factor `ζ(σ+2iγ)` sits near the pole.  The pole's real part is a HELP, not an obstruction: it is
bounded by the same `1/(σ-1)` ceiling uniformly in `γ`.  This is the key analytic input for extending the
zero-free region to small `|t|` (the classical `1-c/log(N(|t|+2))` form / exceptional-zero territory).
EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem hpole_riemannZeta_box :
    ∃ A : ℝ, 0 ≤ A ∧ ∀ s : ℂ, 1 < s.re → s.re ≤ 2 → |s.im| ≤ 1 →
      (-logDeriv riemannZeta s).re ≤ 1 / (s.re - 1) + A := by
  obtain ⟨A, hA0, hA⟩ := exists_logDeriv_LFunctionTrivChar₁_bound_box
  refine ⟨A, hA0, fun s h1 h2 h3 => ?_⟩
  have hs1 : s ≠ 1 := by
    intro h; rw [h] at h1; simp at h1
  have hζ : riemannZeta s ≠ 0 := riemannZeta_ne_zero_of_one_le_re (le_of_lt h1)
  have hdec := logDeriv_LFunctionTrivChar₁_one_eq hs1 hζ
  have hneg : -logDeriv riemannZeta s = 1 / (s - 1) - logDeriv (LFunctionTrivChar₁ 1) s := by
    linear_combination hdec
  have hre1 : (1 / (s - 1)).re ≤ 1 / (s.re - 1) := by
    have hz : 0 < (s - 1).re := by rw [Complex.sub_re, Complex.one_re]; linarith
    have h := re_inv_le_inv_re hz
    rw [Complex.sub_re, Complex.one_re] at h
    rw [one_div]
    exact h
  have habs : -(logDeriv (LFunctionTrivChar₁ 1) s).re ≤ ‖logDeriv (LFunctionTrivChar₁ 1) s‖ := by
    calc -(logDeriv (LFunctionTrivChar₁ 1) s).re
        = (-(logDeriv (LFunctionTrivChar₁ 1) s)).re := by rw [Complex.neg_re]
      _ ≤ ‖-(logDeriv (LFunctionTrivChar₁ 1) s)‖ := Complex.re_le_norm _
      _ = ‖logDeriv (LFunctionTrivChar₁ 1) s‖ := norm_neg _
  calc (-logDeriv riemannZeta s).re
      = (1 / (s - 1)).re - (logDeriv (LFunctionTrivChar₁ 1) s).re := by
        rw [hneg, Complex.sub_re]
    _ ≤ 1 / (s.re - 1) + ‖logDeriv (LFunctionTrivChar₁ 1) s‖ := by linarith [hre1, habs]
    _ ≤ 1 / (s.re - 1) + A := by linarith [hA s (le_of_lt h1) h2 h3]

#print axioms hpole_riemannZeta_box
