import LogDerivTrivChar1OneEq
import SmallTZeroFreeStripT
import Mathlib.Tactic

set_option maxHeartbeats 2000000

open Complex Set Metric DirichletCharacter

/-!
# §1691 — `ζ'/ζ` at SMALL `|t|`: `‖ζ'/ζ(s)‖ ≤ C(T₀) + 1/‖s-1‖` on the zero-free strip

The bounded-height piece of the left-line estimate. Near `s = 1` the only blow-up of `ζ'/ζ` is the
explicit simple pole: by §1570, `ζ'/ζ = logDeriv G - 1/(s-1)` with `G = LFunctionTrivChar₁ 1`
ENTIRE and NONVANISHING on the small-`t` zero-free strip (§ SmallTZeroFreeStripT + Mathlib's
`LFunctionTrivChar₁_apply_one_ne_zero` at the center), so `‖logDeriv G‖` is bounded by a compact-sup
constant there. On the shifted line `Re = 1 - c/log(T+2)` the pole term is `≤ log(T+2)/c` — the
small-`t` stretch contributes `O(log T)`, dominated by the ball regime's `O(log²T)`.
-/

theorem logDeriv_zeta_small_t_bound (T₀ : ℝ) (hT₀ : 0 ≤ T₀) :
    ∃ σ₂ C : ℝ, 1 / 2 ≤ σ₂ ∧ σ₂ < 1 ∧ 0 ≤ C ∧ ∀ s : ℂ, σ₂ < s.re → s.re ≤ 2 →
      |s.im| ≤ T₀ → s ≠ 1 → ‖logDeriv riemannZeta s‖ ≤ C + 1 / ‖s - 1‖ := by
  obtain ⟨σ₀, hσ₀half, hσ₀lt, hstrip⟩ := exists_small_t_zero_free_strip_T T₀ hT₀
  set σ₂ : ℝ := (σ₀ + 1) / 2 with hσ₂def
  have hσ₂gt : σ₀ < σ₂ := by rw [hσ₂def]; linarith
  have hσ₂lt : σ₂ < 1 := by rw [hσ₂def]; linarith
  have hσ₂half : 1 / 2 ≤ σ₂ := by rw [hσ₂def]; linarith
  -- ζ nonvanishing on the closed region K
  set K : Set ℂ := {s : ℂ | σ₂ ≤ s.re ∧ s.re ≤ 2 ∧ |s.im| ≤ T₀} with hKdef
  have hζK : ∀ s ∈ K, s ≠ 1 → riemannZeta s ≠ 0 := by
    intro s hs hs1
    obtain ⟨h1, h2, h3⟩ := hs
    rcases le_total s.re 1 with hre | hre
    · exact hstrip s (lt_of_lt_of_le hσ₂gt h1) hre h3 hs1
    · rcases eq_or_lt_of_le hre with heq | hlt
      · exact hstrip s (lt_of_lt_of_le hσ₂gt h1) heq.ge h3 hs1
      · exact riemannZeta_ne_zero_of_one_lt_re hlt
  -- G is entire and nonvanishing on K
  have hGdiff : Differentiable ℂ (LFunctionTrivChar₁ 1) :=
    differentiable_LFunctionTrivChar₁ 1
  have hGne : ∀ s ∈ K, LFunctionTrivChar₁ 1 s ≠ 0 := by
    intro s hs
    rcases eq_or_ne s 1 with rfl | hs1
    · exact LFunctionTrivChar₁_apply_one_ne_zero 1
    · rw [LFunctionTrivChar₁, Function.update_of_ne hs1]
      apply mul_ne_zero (sub_ne_zero.mpr hs1)
      rw [show LFunctionTrivChar 1 s = riemannZeta s by
        rw [LFunctionTrivChar]
        exact congrFun LFunction_modOne_eq s]
      exact hζK s hs hs1
  -- compactness of K
  have hKcompact : IsCompact K := by
    have hclosed : IsClosed K := by
      rw [hKdef]
      apply IsClosed.inter (isClosed_le continuous_const Complex.continuous_re)
      apply IsClosed.inter (isClosed_le Complex.continuous_re continuous_const)
      have : Continuous fun s : ℂ => |s.im| := Complex.continuous_im.abs
      exact isClosed_le this continuous_const
    have hbounded : Bornology.IsBounded K := by
      apply Bornology.IsBounded.subset (Metric.isBounded_closedBall (x := (0 : ℂ))
        (r := 2 + T₀ + 2))
      intro s hs
      obtain ⟨h1, h2, h3⟩ := hs
      rw [mem_closedBall, dist_zero_right]
      calc ‖s‖ ≤ |s.re| + |s.im| := Complex.norm_le_abs_re_add_abs_im s
        _ ≤ 2 + T₀ + 2 := by
            have hre : |s.re| ≤ 2 := by
              rw [abs_le]
              constructor <;> linarith
            linarith
    exact Metric.isCompact_of_isClosed_isBounded hclosed hbounded
  -- compact sup of ‖logDeriv G‖
  have hGlogcont : ContinuousOn (logDeriv (LFunctionTrivChar₁ 1)) K := by
    intro s hs
    apply ContinuousAt.continuousWithinAt
    rw [logDeriv]
    apply ContinuousAt.div
    · have han : AnalyticAt ℂ (LFunctionTrivChar₁ 1) s := hGdiff.analyticAt s
      exact han.deriv.continuousAt
    · exact (hGdiff s).continuousAt
    · exact hGne s hs
  obtain ⟨C, hC⟩ := hKcompact.exists_bound_of_continuousOn hGlogcont
  refine ⟨σ₂, max C 0, hσ₂half, hσ₂lt, le_max_right _ _, ?_⟩
  intro s hre1 hre2 him hs1
  have hsK : s ∈ K := ⟨hre1.le, hre2, him⟩
  have hζs : riemannZeta s ≠ 0 := hζK s hsK hs1
  have hsplit := logDeriv_LFunctionTrivChar₁_one_eq hs1 hζs
  have hζld : logDeriv riemannZeta s = logDeriv (LFunctionTrivChar₁ 1) s - 1 / (s - 1) := by
    linear_combination -hsplit
  rw [hζld]
  calc ‖logDeriv (LFunctionTrivChar₁ 1) s - 1 / (s - 1)‖
      ≤ ‖logDeriv (LFunctionTrivChar₁ 1) s‖ + ‖1 / (s - 1)‖ := norm_sub_le _ _
    _ ≤ max C 0 + 1 / ‖s - 1‖ := by
        have h1 := hC s hsK
        have h2 : ‖(1 : ℂ) / (s - 1)‖ = 1 / ‖s - 1‖ := by
          rw [norm_div, norm_one]
        rw [h2]
        have := le_max_left C 0
        linarith

#print axioms logDeriv_zeta_small_t_bound
