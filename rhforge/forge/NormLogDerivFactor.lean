import LogDerivFromCodiscreteFactor

set_option maxHeartbeats 1000000

open Complex

/-!
# §1671 — NORM-LANDAU from the factorization: `‖logDeriv ζ z‖ ≤ (Σ|mult|)/δ + ‖logDeriv g z‖`

The quantitative norm form of the pointwise split `logDeriv ζ = Σ d(ρ)/(z-ρ) + logDeriv g`
(§ LogDerivFromCodiscreteFactor): if every zero/pole `ρ` of the divisor keeps distance `≥ δ` from
`z`, the zeros' contribution is at most (total multiplicity)/δ, so
`‖logDeriv ζ z‖ ≤ (Σ_ρ |d ρ|)/δ + ‖logDeriv g z‖`. On the box edge, `δ = c/log T` from the
zero-free region and `‖logDeriv g z‖ ≤ 8M'/R` from §1670 — together the `O(log²T)` bound.
-/

theorem norm_logDeriv_le_of_codiscrete_factor {ζ g : ℂ → ℂ} {d : ℂ → ℤ} {U : Set ℂ} {z : ℂ}
    {δ : ℝ} (hδ : 0 < δ)
    (hU : IsOpen U) (hz : z ∈ U) (hdf : (Function.support d).Finite)
    (h : ζ =ᶠ[Filter.codiscreteWithin U] (∏ᶠ u, (· - u) ^ (d u)) • g)
    (hζc : ContinuousAt (logDeriv ζ) z)
    (hFc : ContinuousAt (logDeriv ((∏ᶠ u, (· - u) ^ (d u)) • g)) z)
    (hdist : ∀ ρ ∈ hdf.toFinset, δ ≤ ‖z - ρ‖)
    (hgz : g z ≠ 0) (hgd : DifferentiableAt ℂ g z) :
    ‖logDeriv ζ z‖ ≤ (∑ ρ ∈ hdf.toFinset, |(d ρ : ℝ)|) / δ + ‖logDeriv g z‖ := by
  have hsub : ∀ ρ ∈ hdf.toFinset, z - ρ ≠ 0 := by
    intro ρ hρ
    have := hdist ρ hρ
    intro h0
    rw [h0, norm_zero] at this
    linarith
  rw [logDeriv_eq_sum_of_codiscrete_factor hU hz hdf h hζc hFc hsub hgz hgd]
  calc ‖(∑ ρ ∈ hdf.toFinset, (d ρ : ℂ) / (z - ρ)) + logDeriv g z‖
      ≤ ‖∑ ρ ∈ hdf.toFinset, (d ρ : ℂ) / (z - ρ)‖ + ‖logDeriv g z‖ := norm_add_le _ _
    _ ≤ (∑ ρ ∈ hdf.toFinset, |(d ρ : ℝ)|) / δ + ‖logDeriv g z‖ := by
        gcongr
        calc ‖∑ ρ ∈ hdf.toFinset, (d ρ : ℂ) / (z - ρ)‖
            ≤ ∑ ρ ∈ hdf.toFinset, ‖(d ρ : ℂ) / (z - ρ)‖ := norm_sum_le _ _
          _ ≤ ∑ ρ ∈ hdf.toFinset, |(d ρ : ℝ)| / δ := by
              apply Finset.sum_le_sum
              intro ρ hρ
              rw [norm_div]
              have hnum : ‖((d ρ : ℤ) : ℂ)‖ = |(d ρ : ℝ)| := by
                rw [show ((d ρ : ℤ) : ℂ) = ((d ρ : ℝ) : ℂ) by push_cast; ring,
                  Complex.norm_real, Real.norm_eq_abs]
              rw [hnum]
              apply div_le_div_of_nonneg_left (abs_nonneg _) hδ (hdist ρ hρ)
          _ = (∑ ρ ∈ hdf.toFinset, |(d ρ : ℝ)|) / δ := by
              rw [Finset.sum_div]

#print axioms norm_logDeriv_le_of_codiscrete_factor
