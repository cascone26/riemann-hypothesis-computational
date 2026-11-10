import ConvexSphere
import ZetaReGeOne

/-!
# S3008 — the combined convexity sphere bound Ms (lane 1)

`convex_sphere_Ms`: ∃ C ≥ 1, for σ₀ > 1, σ₀ − 1 ≤ R ≤ 1, σ₀ + R ≤ 2, |γ| ≥ 3, every z on sphere(σ₀ + iγ, R):
  ‖ζ(z)‖ ≤ C·(2 + |γ|)^{(1 − σ₀ + R)/2}·(log(2 + |γ|) + 6).
(Re z ≤ 1: S3004 convexity; Re z > 1: S3006 log bound, absorbed since the power is ≥ 1.)
⇒ log Ms = ((1−σ₀+R)/2)·log(2+|γ|) + log log + O(1): HALF the sharp N-term exponent — the growth input
(8/R)·log Ms ≈ 4·log|γ| for the next quartic-region rung.
-/

open Complex Metric

namespace ConvexSphereMs

theorem convex_sphere_Ms : ∃ C : ℝ, 1 ≤ C ∧ ∀ (σ₀ γ R : ℝ) (z : ℂ), 1 < σ₀ → σ₀ - 1 ≤ R → R ≤ 1 →
    σ₀ + R ≤ 2 → 3 ≤ |γ| → z ∈ sphere ((σ₀ : ℂ) + (γ : ℂ) * I) R →
    ‖riemannZeta z‖ ≤ C * (2 + |γ|) ^ ((1 - σ₀ + R) / 2) * (Real.log (2 + |γ|) + 6) := by
  obtain ⟨C, hC0, hC⟩ := ConvexSphere.zeta_le_on_sphere_left
  refine ⟨max C 1, le_max_right _ _, fun σ₀ γ R z hσ hR0 hR1 hσR hγ hz => ?_⟩
  have hzc : ‖z - ((σ₀ : ℂ) + (γ : ℂ) * I)‖ = R := by rw [← dist_eq_norm]; exact mem_sphere.mp hz
  have hrel : |z.re - σ₀| ≤ R := by
    have := Complex.abs_re_le_norm (z - ((σ₀ : ℂ) + (γ : ℂ) * I)); rw [hzc] at this; simpa using this
  have him : |z.im - γ| ≤ R := by
    have := Complex.abs_im_le_norm (z - ((σ₀ : ℂ) + (γ : ℂ) * I)); rw [hzc] at this; simpa using this
  rw [abs_le] at hrel
  have hRnn : 0 ≤ R := by linarith
  have hpow1 : 1 ≤ (2 + |γ|) ^ ((1 - σ₀ + R) / 2) :=
    Real.one_le_rpow (by linarith [abs_nonneg γ]) (by linarith)
  have hL6 : 6 ≤ Real.log (2 + |γ|) + 6 := by
    have := Real.log_nonneg (show (1 : ℝ) ≤ 2 + |γ| by linarith [abs_nonneg γ]); linarith
  have hexp : (1 - (σ₀ - R)) / 2 = (1 - σ₀ + R) / 2 := by ring
  rcases le_or_gt z.re 1 with hle | hgt
  · have h := hC σ₀ γ R z hRnn hR1 (by linarith) hz (by linarith [hrel.1]) hle
    rw [hexp] at h
    refine le_trans h ?_
    have : 0 ≤ (2 + |γ|) ^ ((1 - σ₀ + R) / 2) * (Real.log (2 + |γ|) + 6) := by positivity
    have hmC : C ≤ max C 1 := le_max_left _ _
    calc C * (2 + |γ|) ^ ((1 - σ₀ + R) / 2) * (Real.log (2 + |γ|) + 6)
        = C * ((2 + |γ|) ^ ((1 - σ₀ + R) / 2) * (Real.log (2 + |γ|) + 6)) := by ring
      _ ≤ max C 1 * ((2 + |γ|) ^ ((1 - σ₀ + R) / 2) * (Real.log (2 + |γ|) + 6)) :=
          mul_le_mul_of_nonneg_right hmC this
      _ = max C 1 * (2 + |γ|) ^ ((1 - σ₀ + R) / 2) * (Real.log (2 + |γ|) + 6) := by ring
  · have hzim : 2 ≤ |z.im| := by
      have := abs_sub_abs_le_abs_sub γ z.im
      rw [abs_sub_comm] at this; linarith
    have h := ZetaReGeOne.norm_zeta_le_log_right hgt (by linarith [hrel.2]) hzim
    have hzle : |z.im| ≤ 2 + |γ| := by linarith [abs_sub_abs_le_abs_sub z.im γ]
    have hlog : Real.log |z.im| + 6 ≤ Real.log (2 + |γ|) + 6 := by
      have := Real.log_le_log (by linarith) hzle; linarith
    refine le_trans h (le_trans hlog ?_)
    have hm1 : 1 ≤ max C 1 := le_max_right _ _
    have hL0 : 0 ≤ Real.log (2 + |γ|) + 6 := by linarith
    calc Real.log (2 + |γ|) + 6 = 1 * 1 * (Real.log (2 + |γ|) + 6) := by ring
      _ ≤ max C 1 * (2 + |γ|) ^ ((1 - σ₀ + R) / 2) * (Real.log (2 + |γ|) + 6) := by gcongr

end ConvexSphereMs

#print axioms ConvexSphereMs.convex_sphere_Ms
