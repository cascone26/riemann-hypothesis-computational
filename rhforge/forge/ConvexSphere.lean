import ZetaConvexityZeta

/-!
# S3004 — convexity on the pipeline sphere: Re z ≤ 1 part (lane 1: halving log M)

For the Borel–Carathéodory sphere centred at c = σ₀ + iγ (σ₀ > 1, radius R ≤ 1, |γ| ≥ 2), every point z with
0 ≤ Re z ≤ 1 satisfies
  ‖ζ(z)‖ ≤ C·(2 + |γ|)^{(1 − (σ₀ − R))/2}·(log(2 + |γ|) + 6)
(S3002 at height Im z, |Im z − γ| ≤ R ≤ 1 ⇒ 1 + |Im z| ≤ 2 + |γ|, exponent monotone in Re z).
So log of the sphere bound on that part is ((1 − σ₀ + R)/2)·log|γ| + O(log log|γ|) — HALF the sharp N-term exponent
(1 − σ₀ + R) of S2961/S2963.  The Re z > 1 part (log|ζ| = O(log log)) is the companion brick.
-/

open Complex Metric

namespace ConvexSphere

theorem zeta_le_on_sphere_left : ∃ C : ℝ, 0 ≤ C ∧ ∀ (σ₀ γ R : ℝ) (z : ℂ), 0 ≤ R → R ≤ 1 →
    2 ≤ |γ| → z ∈ sphere ((σ₀ : ℂ) + (γ : ℂ) * I) R → 0 ≤ z.re → z.re ≤ 1 →
    ‖riemannZeta z‖ ≤ C * (2 + |γ|) ^ ((1 - (σ₀ - R)) / 2) * (Real.log (2 + |γ|) + 6) := by
  obtain ⟨C, hC0, hC⟩ := ZetaConvexityZeta.zeta_convexity
  refine ⟨C, hC0, fun σ₀ γ R z hR0 hR1 hγ hz hre0 hre1 => ?_⟩
  have hzc : ‖z - ((σ₀ : ℂ) + (γ : ℂ) * I)‖ = R := by rw [← dist_eq_norm]; exact mem_sphere.mp hz
  have him : |z.im - γ| ≤ R := by
    have := Complex.abs_im_le_norm (z - ((σ₀ : ℂ) + (γ : ℂ) * I))
    rw [hzc] at this; simpa using this
  have hrel : |z.re - σ₀| ≤ R := by
    have := Complex.abs_re_le_norm (z - ((σ₀ : ℂ) + (γ : ℂ) * I))
    rw [hzc] at this; simpa using this
  have hzim : 1 ≤ |z.im| := by
    have := abs_sub_abs_le_abs_sub γ z.im
    rw [abs_sub_comm] at this; linarith
  have hzeq : z = ((z.re : ℝ) : ℂ) + ((z.im : ℝ) : ℂ) * I := (Complex.re_add_im z).symm
  have hb := hC z.im z.re hzim hre0 hre1
  rw [← hzeq] at hb
  refine le_trans hb ?_
  have hbase : 1 + |z.im| ≤ 2 + |γ| := by
    linarith [abs_sub_abs_le_abs_sub z.im γ]
  have hexp : (1 - z.re) / 2 ≤ (1 - (σ₀ - R)) / 2 := by
    rw [abs_le] at hrel; linarith [hrel.1]
  have hX1 : 1 ≤ 1 + |z.im| := by linarith [abs_nonneg z.im]
  have hpow : (1 + |z.im|) ^ ((1 - z.re) / 2) ≤ (2 + |γ|) ^ ((1 - (σ₀ - R)) / 2) :=
    le_trans (Real.rpow_le_rpow (by positivity) hbase (by linarith))
      (Real.rpow_le_rpow_of_exponent_le (by linarith [abs_nonneg γ]) hexp)
  have hlog : Real.log (1 + |z.im|) + 6 ≤ Real.log (2 + |γ|) + 6 := by
    have := Real.log_le_log (by positivity) hbase; linarith
  have hl0 : 0 ≤ Real.log (1 + |z.im|) + 6 := by have := Real.log_nonneg hX1; linarith
  have hp0 : 0 ≤ (1 + |z.im|) ^ ((1 - z.re) / 2) := Real.rpow_nonneg (by positivity) _
  calc C * (1 + |z.im|) ^ ((1 - z.re) / 2) * (Real.log (1 + |z.im|) + 6)
      ≤ C * (2 + |γ|) ^ ((1 - (σ₀ - R)) / 2) * (Real.log (2 + |γ|) + 6) := by gcongr

end ConvexSphere

#print axioms ConvexSphere.zeta_le_on_sphere_left
