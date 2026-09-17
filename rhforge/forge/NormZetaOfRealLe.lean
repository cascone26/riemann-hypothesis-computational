import ZetaCrudeBound

/-!
# §1441 — explicit real-axis ζ bound `‖ζ(σ)‖ ≤ 3 + 1/(σ-1)` for `σ > 1`

Just the crude bound §1439 at `N = 1` and real `s = σ`: the partial sum is empty, `1^{anything}=1`,
`‖σ-1‖ = σ-1`, `‖σ‖ = σ`, giving `‖ζ(σ)‖ ≤ 1 + 1 + 1/(σ-1) + σ·(1/σ) = 3 + 1/(σ-1)`.  This explicit upper bound
on `ζ` along the real axis feeds the ζ LOWER bound at `σK+iγ` (via `∑|μ n|/n^σK ≤ ζ(σK)` and §1422), which the BC
scalar `hBexplicit` needs.  EXIT0 / std-axioms.
-/

open Complex

theorem norm_riemannZeta_ofReal_le {σ : ℝ} (hσ : 1 < σ) :
    ‖riemannZeta (σ : ℂ)‖ ≤ 3 + 1 / (σ - 1) := by
  have hσ0 : (0 : ℝ) < (σ : ℂ).re := by rw [Complex.ofReal_re]; linarith
  have hne : (σ : ℂ) ≠ 1 := by rw [Ne, Complex.ofReal_eq_one]; linarith
  have h := norm_riemannZeta_le_crude (N := 1) le_rfl hσ0 hne
  have hnorm1 : ‖(σ : ℂ) - 1‖ = σ - 1 := by
    rw [show (σ : ℂ) - 1 = (((σ - 1 : ℝ)) : ℂ) by push_cast; ring, Complex.norm_real,
      Real.norm_eq_abs, abs_of_pos (by linarith)]
  have hnorms : ‖(σ : ℂ)‖ = σ := by
    rw [Complex.norm_real, Real.norm_eq_abs, abs_of_pos (by linarith)]
  rw [hnorm1, hnorms, Complex.ofReal_re] at h
  simp only [Nat.cast_one, Real.one_rpow] at h
  have hss : σ * (1 / σ) = 1 := by field_simp
  linarith [h, hss]

#print axioms norm_riemannZeta_ofReal_le
