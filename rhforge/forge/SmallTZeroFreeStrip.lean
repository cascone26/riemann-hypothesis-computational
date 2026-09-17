import Mathlib.NumberTheory.LSeries.DirichletContinuation
import Mathlib.NumberTheory.LSeries.Dirichlet
import Mathlib.NumberTheory.LSeries.Nonvanishing

set_option maxHeartbeats 1000000

/-!
# §1573 — SMALL-|t| ZERO-FREE STRIP for ζ:  `∃ σ₀ < 1, no zeros in σ₀ < Re ≤ 1, |Im| ≤ 1`

The first REGION theorem of the small-|t| arc.  Work with the entire function `G = LFunctionTrivChar₁ 1 = (s-1)ζ(s)`
(no pole; `G(1) ≠ 0`).  Its zero set `W` inside the compact rectangle `K₀ = [1/2,1] × [-1,1]` is compact, and `W`
avoids the line `Re = 1` (1-line nonvanishing `riemannZeta_ne_zero_of_one_le_re` + `G(1) ≠ 0`).  If `W = ∅`, take
`σ₀ = 1/2`; else `Re` attains a maximum `σ₀ = Re w* < 1` on `W`, and no ζ-zero can have real part above it.
INEFFECTIVE `σ₀` (compactness, no explicit constant) but GENUINE — the honest classical first step; the effective
`1-c/log` refinement over this strip is the arc's target.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem exists_small_t_zero_free_strip :
    ∃ σ₀ : ℝ, σ₀ < 1 ∧ ∀ s : ℂ, σ₀ < s.re → s.re ≤ 1 → |s.im| ≤ 1 → s ≠ 1 →
      riemannZeta s ≠ 0 := by
  have hG : ∀ s : ℂ, s ≠ 1 → LFunctionTrivChar₁ 1 s = (s - 1) * riemannZeta s := by
    intro s hs
    rw [LFunctionTrivChar₁, Function.update_of_ne hs]
    have h : LFunctionTrivChar 1 s = riemannZeta s := by
      rw [LFunctionTrivChar_eq_mul_riemannZeta hs, Nat.primeFactors_one, Finset.prod_empty, one_mul]
    rw [h]
  have hGcont : Continuous (LFunctionTrivChar₁ 1) :=
    (differentiable_LFunctionTrivChar₁ (n := 1)).continuous
  set K₀ : Set ℂ := Complex.re ⁻¹' Set.Icc (1/2 : ℝ) 1 ∩ Complex.im ⁻¹' Set.Icc (-1 : ℝ) 1 with hK₀def
  have hK₀closed : IsClosed K₀ :=
    (isClosed_Icc.preimage Complex.continuous_re).inter (isClosed_Icc.preimage Complex.continuous_im)
  have hK₀sub : K₀ ⊆ Metric.closedBall (0 : ℂ) 2 := by
    intro s hs
    obtain ⟨hre, him⟩ := hs
    rw [Set.mem_preimage, Set.mem_Icc] at hre him
    rw [Metric.mem_closedBall, dist_zero_right]
    calc ‖s‖ = ‖(s.re : ℂ) + (s.im : ℂ) * I‖ := by rw [Complex.re_add_im]
      _ ≤ ‖((s.re : ℝ) : ℂ)‖ + ‖((s.im : ℝ) : ℂ) * I‖ := norm_add_le _ _
      _ = |s.re| + |s.im| := by
          rw [Complex.norm_real, Real.norm_eq_abs, norm_mul, Complex.norm_I, mul_one,
            Complex.norm_real, Real.norm_eq_abs]
      _ ≤ 2 := by
          have h1 : |s.re| ≤ 1 := abs_le.mpr ⟨by linarith [hre.1], hre.2⟩
          have h2 : |s.im| ≤ 1 := abs_le.mpr ⟨him.1, him.2⟩
          linarith
  have hK₀cpt : IsCompact K₀ :=
    Metric.isCompact_of_isClosed_isBounded hK₀closed (Metric.isBounded_closedBall.subset hK₀sub)
  set W : Set ℂ := K₀ ∩ (LFunctionTrivChar₁ 1) ⁻¹' {0} with hWdef
  have hWcpt : IsCompact W := hK₀cpt.inter_right (isClosed_singleton.preimage hGcont)
  have hWlt1 : ∀ s ∈ W, s.re < 1 := by
    intro s hs
    obtain ⟨hsK, hsz⟩ := hs
    obtain ⟨hre, _⟩ := hsK
    rw [Set.mem_preimage, Set.mem_Icc] at hre
    rw [Set.mem_preimage, Set.mem_singleton_iff] at hsz
    rcases lt_or_eq_of_le hre.2 with h | h
    · exact h
    · exfalso
      by_cases hs1 : s = 1
      · rw [hs1] at hsz
        exact LFunctionTrivChar₁_apply_one_ne_zero 1 hsz
      · rw [hG s hs1] at hsz
        exact (mul_ne_zero (sub_ne_zero.mpr hs1)
          (riemannZeta_ne_zero_of_one_le_re (le_of_eq h.symm))) hsz
  rcases W.eq_empty_or_nonempty with hWe | hWne
  · refine ⟨1/2, by norm_num, fun s h1 h2 h3 hs1 hζ => ?_⟩
    have hsW : s ∈ W := by
      refine ⟨⟨?_, ?_⟩, ?_⟩
      · rw [Set.mem_preimage, Set.mem_Icc]; exact ⟨le_of_lt h1, h2⟩
      · rw [Set.mem_preimage, Set.mem_Icc]; exact abs_le.mp h3
      · rw [Set.mem_preimage, Set.mem_singleton_iff, hG s hs1, hζ, mul_zero]
    rw [hWe] at hsW
    simp at hsW
  · obtain ⟨w, hwW, hwmax⟩ := hWcpt.exists_isMaxOn hWne Complex.continuous_re.continuousOn
    refine ⟨w.re, hWlt1 w hwW, fun s h1 h2 h3 hs1 hζ => ?_⟩
    have hwre : (1/2 : ℝ) ≤ w.re := by
      obtain ⟨hsK, _⟩ := hwW
      obtain ⟨hre, _⟩ := hsK
      rw [Set.mem_preimage, Set.mem_Icc] at hre
      exact hre.1
    have hsW : s ∈ W := by
      refine ⟨⟨?_, ?_⟩, ?_⟩
      · rw [Set.mem_preimage, Set.mem_Icc]; exact ⟨by linarith, h2⟩
      · rw [Set.mem_preimage, Set.mem_Icc]; exact abs_le.mp h3
      · rw [Set.mem_preimage, Set.mem_singleton_iff, hG s hs1, hζ, mul_zero]
    have hle : s.re ≤ w.re := hwmax hsW
    linarith

#print axioms exists_small_t_zero_free_strip
