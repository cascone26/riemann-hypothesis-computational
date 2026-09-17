import Mathlib.NumberTheory.LSeries.DirichletContinuation
import Mathlib.NumberTheory.LSeries.Nonvanishing

set_option maxHeartbeats 1000000

/-!
# §1576 — small-|t| zero-free strip for `LFunction χ` (χ≠1): `∃ σ₀ ∈ [1/2,1), no zeros in σ₀ < Re ≤ 1, |Im| ≤ T`

L-analogue of §1574 — SIMPLER than ζ: `LFunction χ` is ENTIRE for `χ ≠ 1` (no `(s-1)` trick, no `s ≠ 1`
exclusions anywhere).  Zero set `W` in the compact rectangle `[1/2,1] × [-T,T]` is compact; `W` avoids `Re = 1`
by the 1-line nonvanishing `LFunction_ne_zero_of_one_le_re` (with `χ ≠ 1`, valid at every point of the line);
max-Re on `W` gives the ineffective `σ₀ < 1`.  Supplies the strip for the uniform L region.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem exists_small_t_zero_free_strip_LFunction {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) (T : ℝ) (_hT : 0 ≤ T) :
    ∃ σ₀ : ℝ, 1/2 ≤ σ₀ ∧ σ₀ < 1 ∧ ∀ s : ℂ, σ₀ < s.re → s.re ≤ 1 → |s.im| ≤ T →
      LFunction χ s ≠ 0 := by
  have hGcont : Continuous (LFunction χ) := (differentiable_LFunction hχ).continuous
  set K₀ : Set ℂ := Complex.re ⁻¹' Set.Icc (1/2 : ℝ) 1 ∩ Complex.im ⁻¹' Set.Icc (-T : ℝ) T with hK₀def
  have hK₀closed : IsClosed K₀ :=
    (isClosed_Icc.preimage Complex.continuous_re).inter (isClosed_Icc.preimage Complex.continuous_im)
  have hK₀sub : K₀ ⊆ Metric.closedBall (0 : ℂ) (1 + T) := by
    intro s hs
    obtain ⟨hre, him⟩ := hs
    rw [Set.mem_preimage, Set.mem_Icc] at hre him
    rw [Metric.mem_closedBall, dist_zero_right]
    calc ‖s‖ = ‖(s.re : ℂ) + (s.im : ℂ) * I‖ := by rw [Complex.re_add_im]
      _ ≤ ‖((s.re : ℝ) : ℂ)‖ + ‖((s.im : ℝ) : ℂ) * I‖ := norm_add_le _ _
      _ = |s.re| + |s.im| := by
          rw [Complex.norm_real, Real.norm_eq_abs, norm_mul, Complex.norm_I, mul_one,
            Complex.norm_real, Real.norm_eq_abs]
      _ ≤ 1 + T := by
          have h1 : |s.re| ≤ 1 := abs_le.mpr ⟨by linarith [hre.1], hre.2⟩
          have h2 : |s.im| ≤ T := abs_le.mpr ⟨him.1, him.2⟩
          linarith
  have hK₀cpt : IsCompact K₀ :=
    Metric.isCompact_of_isClosed_isBounded hK₀closed (Metric.isBounded_closedBall.subset hK₀sub)
  set W : Set ℂ := K₀ ∩ (LFunction χ) ⁻¹' {0} with hWdef
  have hWcpt : IsCompact W := hK₀cpt.inter_right (isClosed_singleton.preimage hGcont)
  have hWlt1 : ∀ s ∈ W, s.re < 1 := by
    intro s hs
    obtain ⟨hsK, hsz⟩ := hs
    obtain ⟨hre, _⟩ := hsK
    rw [Set.mem_preimage, Set.mem_Icc] at hre
    rw [Set.mem_preimage, Set.mem_singleton_iff] at hsz
    rcases lt_or_eq_of_le hre.2 with h | h
    · exact h
    · exact absurd hsz (LFunction_ne_zero_of_one_le_re χ (Or.inl hχ) (le_of_eq h.symm))
  rcases W.eq_empty_or_nonempty with hWe | hWne
  · refine ⟨1/2, le_refl _, by norm_num, fun s h1 h2 h3 hζ => ?_⟩
    have hsW : s ∈ W := by
      refine ⟨⟨?_, ?_⟩, ?_⟩
      · rw [Set.mem_preimage, Set.mem_Icc]; exact ⟨le_of_lt h1, h2⟩
      · rw [Set.mem_preimage, Set.mem_Icc]; exact abs_le.mp h3
      · rw [Set.mem_preimage, Set.mem_singleton_iff]; exact hζ
    rw [hWe] at hsW
    simp at hsW
  · obtain ⟨w, hwW, hwmax⟩ := hWcpt.exists_isMaxOn hWne Complex.continuous_re.continuousOn
    have hwre : (1/2 : ℝ) ≤ w.re := by
      obtain ⟨hsK, _⟩ := hwW
      obtain ⟨hre, _⟩ := hsK
      rw [Set.mem_preimage, Set.mem_Icc] at hre
      exact hre.1
    refine ⟨w.re, hwre, hWlt1 w hwW, fun s h1 h2 h3 hζ => ?_⟩
    have hsW : s ∈ W := by
      refine ⟨⟨?_, ?_⟩, ?_⟩
      · rw [Set.mem_preimage, Set.mem_Icc]; exact ⟨by linarith, h2⟩
      · rw [Set.mem_preimage, Set.mem_Icc]; exact abs_le.mp h3
      · rw [Set.mem_preimage, Set.mem_singleton_iff]; exact hζ
    have hle : s.re ≤ w.re := hwmax hsW
    linarith

#print axioms exists_small_t_zero_free_strip_LFunction
