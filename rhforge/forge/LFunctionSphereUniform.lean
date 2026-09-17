import NormLFunctionCrudeRePos

/-!
# §1534 — explicit `‖LFunction χ‖` bound uniform over a sphere `sphere c R`

L-analogue of ζ §1440 (`norm_riemannZeta_le_sphere_uniform`), but SIMPLER: `LFunction χ` (χ≠1) is entire, so there is
NO pole term — the crude bound §1533 has only three pieces `N + (N+1)^{-σ}·M + ‖s‖·N^{-σ}/σ·M`.  Each is bounded by
its worst case over the sphere via `z.re ≥ c.re-R` (so `σ` is bounded below), `‖z‖ ≤ ‖c‖+R` (from `‖z-c‖=R` +
triangle).  Since `N,N+1 ≥ 1`, `t ↦ (·)^{-σ}` decreases in `σ`, so the smallest `σ = c.re-R` gives the largest
value.  Yields a SINGLE explicit constant `Ms` bounding `‖LFunction χ‖` on the whole sphere — the supplier for the
`hLsphere` hypothesis of §1500 once `c = σK+iγ`, `R`, `N` are fixed.  Choosing `N = ⌈|γ|⌉` makes `Ms` polynomial in
`|γ|`, so `log Ms ≍ log|γ|`.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric

theorem norm_LFunction_le_sphere_uniform {M : ℕ} [NeZero M] {χ : DirichletCharacter ℂ M} (hχ : χ ≠ 1)
    {c : ℂ} {R : ℝ} {N : ℕ} (hN : 1 ≤ N) (hR : 0 < R) (hcR : 0 < c.re - R) {z : ℂ}
    (hz : z ∈ sphere c R) :
    ‖LFunction χ z‖ ≤ (N : ℝ)
      + (((N + 1 : ℕ) : ℝ) ^ (-(c.re - R)) * M
        + (‖c‖ + R) * (N : ℝ) ^ (-(c.re - R)) / (c.re - R) * M) := by
  have hNR : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast hN
  have hN1 : (1 : ℝ) ≤ ((N + 1 : ℕ) : ℝ) := by exact_mod_cast Nat.le_add_left 1 N
  have hzc : ‖z - c‖ = R := by rw [← dist_eq_norm]; exact (mem_sphere.mp hz)
  have hzre : c.re - R ≤ z.re := by
    have h : |z.re - c.re| ≤ R := by
      rw [← hzc]; simpa [Complex.sub_re] using Complex.abs_re_le_norm (z - c)
    rw [abs_le] at h; linarith [h.1]
  have hzn : ‖z‖ ≤ ‖c‖ + R := by
    have : z = c + (z - c) := by ring
    rw [this]; exact le_trans (norm_add_le _ _) (by rw [hzc])
  have hzre0 : 0 < z.re := lt_of_lt_of_le hcR hzre
  have hcrude := norm_LFunction_le_crude_re_pos hχ hzre0 hN
  calc ‖LFunction χ z‖
      ≤ (N : ℝ) + (((N + 1 : ℕ) : ℝ) ^ (-z.re) * M + ‖z‖ * (N : ℝ) ^ (-z.re) / z.re * M) := hcrude
    _ ≤ (N : ℝ) + (((N + 1 : ℕ) : ℝ) ^ (-(c.re - R)) * M
          + (‖c‖ + R) * (N : ℝ) ^ (-(c.re - R)) / (c.re - R) * M) := by
        have hMnn : (0 : ℝ) ≤ (M : ℝ) := by positivity
        gcongr

#print axioms norm_LFunction_le_sphere_uniform
