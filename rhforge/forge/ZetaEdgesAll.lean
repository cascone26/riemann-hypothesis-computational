import ZetaEdgesUniform

/-!
# S2993 — uniform convexity edges for ALL real t (lane 1)

Compactness closes |t| ≤ 2 (E = (s−1)ζ(s) is entire, S2979), S2991 covers |t| ≥ 2:
* `exists_bound_small` — ∃ M ≥ 0, |t| ≤ 2 ⇒ ‖E(it)‖ ≤ M ∧ ‖E(1+it)‖ ≤ M.
* `left_edge_all`  — ∃ K, ∀ t, ‖E(it)‖   ≤ K·(1+|t|)^{3/2}·(log(1+|t|) + 5).
* `right_edge_all` — ∃ K, ∀ t, ‖E(1+it)‖ ≤ K·(1+|t|)·(log(1+|t|) + 5).
Exactly the shape consumed by S2986/S2988 (Gaussian localisation with exponent p = 3/2, 1) ⇒ explicit a, b for S2984.
-/

open Complex Set

namespace ZetaEdgesAll

theorem exists_bound_small : ∃ M : ℝ, 0 ≤ M ∧ ∀ t : ℝ, |t| ≤ 2 →
    ‖ZetaEntire.zetaE ((t : ℂ) * I)‖ ≤ M ∧ ‖ZetaEntire.zetaE (1 + (t : ℂ) * I)‖ ≤ M := by
  have hc1 : Continuous fun t : ℝ => ZetaEntire.zetaE ((t : ℂ) * I) :=
    ZetaEntire.zetaE_differentiable.continuous.comp (by fun_prop)
  have hc2 : Continuous fun t : ℝ => ZetaEntire.zetaE (1 + (t : ℂ) * I) :=
    ZetaEntire.zetaE_differentiable.continuous.comp (by fun_prop)
  obtain ⟨M1, hM1⟩ := (isCompact_Icc (a := (-2 : ℝ)) (b := 2)).exists_bound_of_continuousOn
    hc1.continuousOn
  obtain ⟨M2, hM2⟩ := (isCompact_Icc (a := (-2 : ℝ)) (b := 2)).exists_bound_of_continuousOn
    hc2.continuousOn
  refine ⟨max (max M1 M2) 0, le_max_right _ _, fun t ht => ?_⟩
  have hmem : t ∈ Icc (-2 : ℝ) 2 := by rw [abs_le] at ht; exact ⟨ht.1, ht.2⟩
  exact ⟨le_trans (hM1 t hmem) (le_trans (le_max_left _ _) (le_max_left _ _)),
    le_trans (hM2 t hmem) (le_trans (le_max_right _ _) (le_max_left _ _))⟩

theorem left_edge_all : ∃ K : ℝ, ∀ t : ℝ,
    ‖ZetaEntire.zetaE ((t : ℂ) * I)‖ ≤ K * (1 + |t|) ^ (3 / 2 : ℝ) * (Real.log (1 + |t|) + 5) := by
  obtain ⟨M, hM0, hM⟩ := exists_bound_small
  refine ⟨max M 1, fun t => ?_⟩
  have hK1 : (1 : ℝ) ≤ max M 1 := le_max_right _ _
  have hp1 : (1 : ℝ) ≤ (1 + |t|) ^ (3 / 2 : ℝ) :=
    Real.one_le_rpow (by linarith [abs_nonneg t]) (by norm_num)
  have hlog0 : 0 ≤ Real.log (1 + |t|) := Real.log_nonneg (by linarith [abs_nonneg t])
  rcases le_or_gt |t| 2 with hs | hl
  · have h := (hM t hs).1
    calc ‖ZetaEntire.zetaE ((t : ℂ) * I)‖ ≤ M := h
      _ ≤ max M 1 * 1 * 1 := by simp
      _ ≤ max M 1 * (1 + |t|) ^ (3 / 2 : ℝ) * (Real.log (1 + |t|) + 5) := by
          gcongr; linarith
  · have hne : ((t : ℂ) * I) ≠ 1 := by
      intro h; have := congrArg Complex.re h; simp at this
    rw [ZetaEntire.zetaE_eq hne, norm_mul]
    have h1 : ‖(t : ℂ) * I - 1‖ ≤ 1 + |t| := by
      calc ‖(t : ℂ) * I - 1‖ ≤ ‖(t : ℂ) * I‖ + ‖(1 : ℂ)‖ := norm_sub_le _ _
        _ = 1 + |t| := by
          rw [norm_one, norm_mul, Complex.norm_I, mul_one, Complex.norm_real, Real.norm_eq_abs]
          ring
    have h2 := ZetaEdgesUniform.norm_zeta_mul_I_le_abs hl.le
    have hsq : Real.sqrt (|t| / (2 * Real.pi)) ≤ (1 + |t|) ^ (1 / 2 : ℝ) := by
      rw [← Real.sqrt_eq_rpow]
      apply Real.sqrt_le_sqrt
      have : |t| / (2 * Real.pi) ≤ |t| := by
        rw [div_le_iff₀ (by positivity)]
        nlinarith [Real.pi_gt_three, abs_nonneg t]
      linarith
    have hlogle : Real.log |t| ≤ Real.log (1 + |t|) :=
      Real.log_le_log (by linarith) (by linarith)
    have hpow : (1 + |t|) * (1 + |t|) ^ (1 / 2 : ℝ) = (1 + |t|) ^ (3 / 2 : ℝ) := by
      rw [show (3 / 2 : ℝ) = 1 + 1 / 2 by norm_num, Real.rpow_add (by positivity), Real.rpow_one]
    have hlog5 : 0 ≤ Real.log |t| + 5 := by
      have := Real.log_nonneg (show (1 : ℝ) ≤ |t| by linarith); linarith
    calc ‖(t : ℂ) * I - 1‖ * ‖riemannZeta ((t : ℂ) * I)‖
        ≤ (1 + |t|) * (Real.sqrt (|t| / (2 * Real.pi)) * (Real.log |t| + 5)) :=
          mul_le_mul h1 h2 (norm_nonneg _) (by positivity)
      _ ≤ (1 + |t|) * ((1 + |t|) ^ (1 / 2 : ℝ) * (Real.log (1 + |t|) + 5)) := by
          gcongr
      _ = (1 + |t|) ^ (3 / 2 : ℝ) * (Real.log (1 + |t|) + 5) := by rw [← hpow]; ring
      _ ≤ max M 1 * (1 + |t|) ^ (3 / 2 : ℝ) * (Real.log (1 + |t|) + 5) := by
          have hA : 0 ≤ (1 + |t|) ^ (3 / 2 : ℝ) * (Real.log (1 + |t|) + 5) := by positivity
          nlinarith

theorem right_edge_all : ∃ K : ℝ, ∀ t : ℝ,
    ‖ZetaEntire.zetaE (1 + (t : ℂ) * I)‖ ≤ K * (1 + |t|) * (Real.log (1 + |t|) + 5) := by
  obtain ⟨M, hM0, hM⟩ := exists_bound_small
  refine ⟨max M 1, fun t => ?_⟩
  have hK1 : (1 : ℝ) ≤ max M 1 := le_max_right _ _
  have hlog0 : 0 ≤ Real.log (1 + |t|) := Real.log_nonneg (by linarith [abs_nonneg t])
  rcases le_or_gt |t| 2 with hs | hl
  · have h := (hM t hs).2
    calc ‖ZetaEntire.zetaE (1 + (t : ℂ) * I)‖ ≤ M := h
      _ ≤ max M 1 * 1 * 1 := by simp
      _ ≤ max M 1 * (1 + |t|) * (Real.log (1 + |t|) + 5) := by
          gcongr
          · linarith [abs_nonneg t]
          · linarith
  · have hne : (1 + (t : ℂ) * I) ≠ 1 := by
      intro h; have := congrArg Complex.im h; simp at this
      rw [this] at hl; norm_num at hl
    rw [ZetaEntire.zetaE_eq hne, norm_mul,
      show (1 + (t : ℂ) * I) - 1 = (t : ℂ) * I by ring, norm_mul, Complex.norm_I, mul_one,
      Complex.norm_real, Real.norm_eq_abs]
    have h2 := ZetaEdgesUniform.norm_zeta_one_add_mul_I_le_abs hl.le
    have hlogle : Real.log |t| ≤ Real.log (1 + |t|) :=
      Real.log_le_log (by linarith) (by linarith)
    calc |t| * ‖riemannZeta (1 + (t : ℂ) * I)‖ ≤ (1 + |t|) * (Real.log (1 + |t|) + 5) := by
          apply mul_le_mul (by linarith) (by linarith) (norm_nonneg _) (by linarith)
      _ ≤ max M 1 * (1 + |t|) * (Real.log (1 + |t|) + 5) := by
          have hA : 0 ≤ (1 + |t|) * (Real.log (1 + |t|) + 5) := by positivity
          nlinarith

end ZetaEdgesAll

#print axioms ZetaEdgesAll.exists_bound_small
#print axioms ZetaEdgesAll.left_edge_all
#print axioms ZetaEdgesAll.right_edge_all
