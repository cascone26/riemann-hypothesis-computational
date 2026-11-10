import ZetaThreeLinesHB
import ZetaEdgesAll
import GaussWeightLog

/-!
# S3000 — THE CONVEXITY BOUND for ζ, fully formal (lane 1)

`zetaE_convexity`: ∃ C, ∀ T, ∀ σ ∈ [0,1],
  ‖(s−1)ζ(s)‖ ≤ C·(1+|T|)^{(3−σ)/2}·(log(1+|T|) + 6)   at s = σ + iT,
i.e. |ζ(σ+iT)| ≪ |T|^{(1−σ)/2}·log|T| — the classical convexity (Lindelöf-type) bound.
Chain: E entire (S2979) · exact |ζ(it)| via functional equation (S2975) · edges for all t (S2977/S2982/S2991/S2993) ·
Gaussian localiser + Hadamard three-lines (S2984) · exact-exponent Gaussian moments with log (S2986/S2988) ·
strip growth through 1/Γ reflection (S2995/S2997) · hB (S2999).
-/

open Complex

namespace ZetaConvexity

theorem rpow_mix_le {A B : ℝ} (hA : 0 ≤ A) (hB : 0 ≤ B) {σ : ℝ} (h0 : 0 ≤ σ) (h1 : σ ≤ 1) :
    A ^ (1 - σ) * B ^ σ ≤ max A B := by
  have hM : 0 ≤ max A B := le_trans hA (le_max_left _ _)
  have e1 : A ^ (1 - σ) ≤ (max A B) ^ (1 - σ) :=
    Real.rpow_le_rpow hA (le_max_left _ _) (by linarith)
  have e2 : B ^ σ ≤ (max A B) ^ σ := Real.rpow_le_rpow hB (le_max_right _ _) h0
  calc A ^ (1 - σ) * B ^ σ ≤ (max A B) ^ (1 - σ) * (max A B) ^ σ :=
        mul_le_mul e1 e2 (Real.rpow_nonneg hB _) (Real.rpow_nonneg hM _)
    _ = max A B := by
        rcases eq_or_lt_of_le hM with h | h
        · rw [← h]
          rcases eq_or_lt_of_le h0 with hs | hs
          · subst hs; simp
          · rw [Real.zero_rpow hs.ne', mul_zero]
        · rw [← Real.rpow_add h]; simp

theorem zetaE_convexity : ∃ C : ℝ, 0 ≤ C ∧ ∀ T σ : ℝ, 0 ≤ σ → σ ≤ 1 →
    ‖ZetaEntire.zetaE ((σ : ℂ) + (T : ℂ) * I)‖
      ≤ C * (1 + |T|) ^ ((3 - σ) / 2) * (Real.log (1 + |T|) + 6) := by
  obtain ⟨K1, hK1⟩ := ZetaEdgesAll.left_edge_all
  obtain ⟨K2, hK2⟩ := ZetaEdgesAll.right_edge_all
  set k1 := max K1 0 with hk1
  set k2 := max K2 0 with hk2
  have hk10 : 0 ≤ k1 := le_max_right _ _
  have hk20 : 0 ≤ k2 := le_max_right _ _
  set A := k1 * Real.exp (25 / 16) with hA
  set B := Real.exp 1 * k2 * Real.exp 1 with hB
  have hA0 : 0 ≤ A := by positivity
  have hB0 : 0 ≤ B := by positivity
  refine ⟨max A B, le_trans hA0 (le_max_left _ _), fun T σ h0 h1 => ?_⟩
  set X := 1 + |T| with hX
  set Lg := Real.log (1 + |T|) + 6 with hLg
  have hX1 : 1 ≤ X := by rw [hX]; linarith [abs_nonneg T]
  have hX0 : 0 < X := by linarith
  have hLg0 : 0 ≤ Lg := by rw [hLg]; have := Real.log_nonneg hX1; linarith
  -- localised left edge
  have ha : ∀ t : ℝ, ‖ZetaEntire.zetaE ((t : ℂ) * I)‖ * Real.exp (-(t - T) ^ 2)
      ≤ A * X ^ (3 / 2 : ℝ) * Lg := by
    intro t
    have hE := hK1 t
    have hE' : ‖ZetaEntire.zetaE ((t : ℂ) * I)‖
        ≤ k1 * (1 + |t|) ^ (3 / 2 : ℝ) * (Real.log (1 + |t|) + 5) := by
      refine le_trans hE ?_
      have : 0 ≤ (1 + |t|) ^ (3 / 2 : ℝ) * (Real.log (1 + |t|) + 5) := by
        have := Real.log_nonneg (show (1 : ℝ) ≤ 1 + |t| by linarith [abs_nonneg t]); positivity
      nlinarith [le_max_left K1 0]
    have hw := GaussWeightLog.gauss_poly_log_weight t T (p := 3 / 2) (c := 5) (by norm_num)
      (by norm_num)
    have hg0 : 0 ≤ Real.exp (-(t - T) ^ 2) := (Real.exp_pos _).le
    calc ‖ZetaEntire.zetaE ((t : ℂ) * I)‖ * Real.exp (-(t - T) ^ 2)
        ≤ k1 * (1 + |t|) ^ (3 / 2 : ℝ) * (Real.log (1 + |t|) + 5) * Real.exp (-(t - T) ^ 2) :=
          mul_le_mul_of_nonneg_right hE' hg0
      _ = k1 * ((1 + |t|) ^ (3 / 2 : ℝ) * (Real.log (1 + |t|) + 5) * Real.exp (-(t - T) ^ 2)) := by
          ring
      _ ≤ k1 * ((1 + |T|) ^ (3 / 2 : ℝ) * (Real.log (1 + |T|) + 5 + 1)
            * Real.exp ((3 / 2 + 1) ^ 2 / 4)) := mul_le_mul_of_nonneg_left hw hk10
      _ = A * X ^ (3 / 2 : ℝ) * Lg := by
          rw [hA, hX, hLg]; norm_num; ring
  -- localised right edge
  have hb : ∀ t : ℝ, ‖ZetaEntire.zetaE (1 + (t : ℂ) * I)‖ * Real.exp (1 - (t - T) ^ 2)
      ≤ B * X * Lg := by
    intro t
    have hE := hK2 t
    have hE' : ‖ZetaEntire.zetaE (1 + (t : ℂ) * I)‖
        ≤ k2 * (1 + |t|) ^ (1 : ℝ) * (Real.log (1 + |t|) + 5) := by
      rw [Real.rpow_one]
      refine le_trans hE ?_
      have : 0 ≤ (1 + |t|) * (Real.log (1 + |t|) + 5) := by
        have := Real.log_nonneg (show (1 : ℝ) ≤ 1 + |t| by linarith [abs_nonneg t]); positivity
      nlinarith [le_max_left K2 0]
    have hw := GaussWeightLog.gauss_poly_log_weight t T (p := 1) (c := 5) (by norm_num) (by norm_num)
    have hg0 : 0 ≤ Real.exp (-(t - T) ^ 2) := (Real.exp_pos _).le
    have hsplit : Real.exp (1 - (t - T) ^ 2) = Real.exp 1 * Real.exp (-(t - T) ^ 2) := by
      rw [← Real.exp_add]; ring_nf
    rw [hsplit]
    calc ‖ZetaEntire.zetaE (1 + (t : ℂ) * I)‖ * (Real.exp 1 * Real.exp (-(t - T) ^ 2))
        ≤ k2 * (1 + |t|) ^ (1 : ℝ) * (Real.log (1 + |t|) + 5)
            * (Real.exp 1 * Real.exp (-(t - T) ^ 2)) :=
          mul_le_mul_of_nonneg_right hE' (by positivity)
      _ = Real.exp 1 * k2 * ((1 + |t|) ^ (1 : ℝ) * (Real.log (1 + |t|) + 5)
            * Real.exp (-(t - T) ^ 2)) := by ring
      _ ≤ Real.exp 1 * k2 * ((1 + |T|) ^ (1 : ℝ) * (Real.log (1 + |T|) + 5 + 1)
            * Real.exp ((1 + 1) ^ 2 / 4)) := mul_le_mul_of_nonneg_left hw (by positivity)
      _ = B * X * Lg := by
          rw [hB, hX, hLg, Real.rpow_one]; norm_num; ring
  have key := ZetaThreeLines.zetaE_three_lines T (A * X ^ (3 / 2 : ℝ) * Lg) (B * X * Lg) σ h0 h1
    (ZetaThreeLinesHB.zetaE_gauss_bddAbove T) ha hb
  have hσ2 : 1 ≤ Real.exp (σ ^ 2) := Real.one_le_exp (sq_nonneg σ)
  have hEn : 0 ≤ ‖ZetaEntire.zetaE ((σ : ℂ) + (T : ℂ) * I)‖ := norm_nonneg _
  have hstep : ‖ZetaEntire.zetaE ((σ : ℂ) + (T : ℂ) * I)‖
      ≤ (A * X ^ (3 / 2 : ℝ) * Lg) ^ (1 - σ) * (B * X * Lg) ^ σ := by
    nlinarith
  have hXa : 0 ≤ X ^ (3 / 2 : ℝ) := Real.rpow_nonneg hX0.le _
  have hL : Lg ^ (1 - σ) * Lg ^ σ = Lg := by
    rcases eq_or_lt_of_le hLg0 with h | h
    · rw [← h]
      rcases eq_or_lt_of_le h0 with hs | hs
      · subst hs; simp
      · rw [Real.zero_rpow hs.ne']; simp
    · rw [← Real.rpow_add h]; simp
  have hXp : X ^ (3 / 2 * (1 - σ)) * X ^ σ = X ^ ((3 - σ) / 2) := by
    rw [← Real.rpow_add hX0]; congr 1; ring
  have hmix := rpow_mix_le hA0 hB0 h0 h1
  have hXr : 0 ≤ X ^ ((3 - σ) / 2) := Real.rpow_nonneg hX0.le _
  have heq : (A * X ^ (3 / 2 : ℝ) * Lg) ^ (1 - σ) * (B * X * Lg) ^ σ
      = A ^ (1 - σ) * B ^ σ * X ^ ((3 - σ) / 2) * Lg := by
    rw [Real.mul_rpow (by positivity) hLg0, Real.mul_rpow hA0 hXa, Real.mul_rpow (by positivity) hLg0,
      Real.mul_rpow hB0 hX0.le, ← Real.rpow_mul hX0.le]
    calc A ^ (1 - σ) * X ^ (3 / 2 * (1 - σ)) * Lg ^ (1 - σ) * (B ^ σ * X ^ σ * Lg ^ σ)
        = A ^ (1 - σ) * B ^ σ * (X ^ (3 / 2 * (1 - σ)) * X ^ σ) * (Lg ^ (1 - σ) * Lg ^ σ) := by ring
      _ = A ^ (1 - σ) * B ^ σ * X ^ ((3 - σ) / 2) * Lg := by rw [hXp, hL]
  calc ‖ZetaEntire.zetaE ((σ : ℂ) + (T : ℂ) * I)‖
      ≤ (A * X ^ (3 / 2 : ℝ) * Lg) ^ (1 - σ) * (B * X * Lg) ^ σ := hstep
    _ = A ^ (1 - σ) * B ^ σ * X ^ ((3 - σ) / 2) * Lg := heq
    _ ≤ max A B * X ^ ((3 - σ) / 2) * Lg := by gcongr

end ZetaConvexity

#print axioms ZetaConvexity.rpow_mix_le
#print axioms ZetaConvexity.zetaE_convexity
