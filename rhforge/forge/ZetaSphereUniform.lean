import ZetaCrudeBound

/-!
# §1440 — explicit ζ bound uniform over a sphere `sphere c R`

Bounds every term of the crude ζ bound §1439 by its worst case over the sphere: `z.re ≥ c.re-R`,
`‖z-1‖ ≥ ‖c-1‖-R`, `‖z‖ ≤ ‖c‖+R` (from `‖z-c‖ = R` + triangle inequalities).  Gives a SINGLE explicit constant
`Ms` bounding `‖ζ‖` on the whole sphere — the supplier for §1437's `hζsphere` once `c = σK+iγ`, `R`, `N` are fixed.
EXIT0 / std-axioms.
-/

open Complex Metric

theorem norm_riemannZeta_le_sphere_uniform {c : ℂ} {R : ℝ} {N : ℕ} (hN : 1 ≤ N) (hR : 0 < R)
    (hcR : 0 < c.re - R) (hc1 : R < ‖c - 1‖) {z : ℂ} (hz : z ∈ sphere c R) :
    ‖riemannZeta z‖ ≤ (N : ℝ) + 1
      + (N : ℝ) ^ (1 - (c.re - R)) / (‖c - 1‖ - R)
      + (‖c‖ + R) * ((N : ℝ) ^ (-(c.re - R)) / (c.re - R)) := by
  have hNR : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast hN
  have hzc : ‖z - c‖ = R := by rw [← dist_eq_norm]; exact (mem_sphere.mp hz)
  have hzre : c.re - R ≤ z.re := by
    have h : |z.re - c.re| ≤ R := by
      rw [← hzc]; simpa [Complex.sub_re] using Complex.abs_re_le_norm (z - c)
    rw [abs_le] at h; linarith [h.1]
  have hz1 : ‖c - 1‖ - R ≤ ‖z - 1‖ := by
    have htri : ‖c - 1‖ ≤ ‖c - z‖ + ‖z - 1‖ := by
      have : c - 1 = (c - z) + (z - 1) := by ring
      rw [this]; exact norm_add_le _ _
    have hcz : ‖c - z‖ = R := by rw [norm_sub_rev]; exact hzc
    linarith
  have hzn : ‖z‖ ≤ ‖c‖ + R := by
    have : z = c + (z - c) := by ring
    rw [this]; exact le_trans (norm_add_le _ _) (by rw [hzc])
  have hzre0 : 0 < z.re := lt_of_lt_of_le hcR hzre
  have hz1pos : 0 < ‖c - 1‖ - R := by linarith
  have hzne1 : z ≠ 1 := by
    intro h; rw [h] at hz1; simp at hz1; linarith
  have hcrude := norm_riemannZeta_le_crude hN hzre0 hzne1
  calc ‖riemannZeta z‖
      ≤ (N : ℝ) + 1 + (N : ℝ) ^ (1 - z.re) / ‖z - 1‖ + ‖z‖ * ((N : ℝ) ^ (-z.re) / z.re) := hcrude
    _ ≤ (N : ℝ) + 1 + (N : ℝ) ^ (1 - (c.re - R)) / (‖c - 1‖ - R)
          + (‖c‖ + R) * ((N : ℝ) ^ (-(c.re - R)) / (c.re - R)) := by
        gcongr

#print axioms norm_riemannZeta_le_sphere_uniform
