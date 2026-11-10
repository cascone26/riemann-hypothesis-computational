import ZetaConvexity

/-!
# S3002 — convexity bound for ζ itself (lane 1: input for the sphere bound of the zero-free pipeline)

`zeta_convexity`: ∃ C ≥ 0, ∀ T with |T| ≥ 1, ∀ σ ∈ [0,1]:
  ‖ζ(σ + iT)‖ ≤ C·(1+|T|)^{(1−σ)/2}·(log(1+|T|) + 6).
From S3000 (E = (s−1)ζ(s)) dividing by ‖s − 1‖ ≥ |T| ≥ (1+|T|)/2 and (1+|T|)^{(3−σ)/2} = (1+|T|)·(1+|T|)^{(1−σ)/2}.
Use in the pipeline: on sphere(σ₀+iγ, R) the points with Re ≤ 1 have log|ζ| ≤ ((1−Re)/2)·log|γ| + O(log log|γ|) ≤ (R/2)·log|γ| + …,
halving log M vs the sharp N-term bound (S2963: (1−σ+R)·log) ⇒ leading growth coefficient 8 → 4 (Report until wired).
-/

open Complex

namespace ZetaConvexityZeta

theorem zeta_convexity : ∃ C : ℝ, 0 ≤ C ∧ ∀ T σ : ℝ, 1 ≤ |T| → 0 ≤ σ → σ ≤ 1 →
    ‖riemannZeta ((σ : ℂ) + (T : ℂ) * I)‖
      ≤ C * (1 + |T|) ^ ((1 - σ) / 2) * (Real.log (1 + |T|) + 6) := by
  obtain ⟨C, hC0, hC⟩ := ZetaConvexity.zetaE_convexity
  refine ⟨2 * C, by positivity, fun T σ hT h0 h1 => ?_⟩
  set s : ℂ := (σ : ℂ) + (T : ℂ) * I with hs
  have hsne : s ≠ 1 := by
    intro h; have := congrArg Complex.im h; simp [hs] at this; rw [this] at hT; norm_num at hT
  have hE := hC T σ h0 h1
  rw [ZetaEntire.zetaE_eq hsne, norm_mul] at hE
  have hT1 : |T| ≤ ‖s - 1‖ := by
    have := Complex.abs_im_le_norm (s - 1); simpa [hs] using this
  have hX : 0 < 1 + |T| := by positivity
  have hLg : 0 ≤ Real.log (1 + |T|) + 6 := by
    have := Real.log_nonneg (show (1 : ℝ) ≤ 1 + |T| by linarith [abs_nonneg T]); linarith
  have hsplit : (1 + |T|) ^ ((3 - σ) / 2) = (1 + |T|) * (1 + |T|) ^ ((1 - σ) / 2) := by
    rw [show (3 - σ) / 2 = 1 + (1 - σ) / 2 by ring, Real.rpow_add hX, Real.rpow_one]
  rw [hsplit] at hE
  set Y := (1 + |T|) ^ ((1 - σ) / 2) with hY
  have hY0 : 0 ≤ Y := Real.rpow_nonneg hX.le _
  have hz0 : 0 ≤ ‖riemannZeta s‖ := norm_nonneg _
  -- ‖s−1‖·‖ζ‖ ≤ C(1+|T|)Y·Lg and ‖s−1‖ ≥ |T| ≥ (1+|T|)/2
  have hhalf : (1 + |T|) / 2 ≤ ‖s - 1‖ := by linarith
  have hprod : (1 + |T|) / 2 * ‖riemannZeta s‖ ≤ C * ((1 + |T|) * Y) * (Real.log (1 + |T|) + 6) :=
    le_trans (mul_le_mul_of_nonneg_right hhalf hz0) hE
  have hpos : 0 < (1 + |T|) / 2 := by positivity
  have key : (1 + |T|) / 2 * ‖riemannZeta s‖ ≤ (1 + |T|) / 2 * (2 * C * Y * (Real.log (1 + |T|) + 6)) := by
    calc (1 + |T|) / 2 * ‖riemannZeta s‖ ≤ C * ((1 + |T|) * Y) * (Real.log (1 + |T|) + 6) := hprod
      _ = (1 + |T|) / 2 * (2 * C * Y * (Real.log (1 + |T|) + 6)) := by ring
  exact le_of_mul_le_mul_left key hpos

end ZetaConvexityZeta

#print axioms ZetaConvexityZeta.zeta_convexity
