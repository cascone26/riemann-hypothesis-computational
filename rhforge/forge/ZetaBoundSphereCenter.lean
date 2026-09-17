import ZetaBoundOnBall
import SphereSubsetBall

/-!
# §1398 — ζ-bound on the sphere centred at the test point `σ+iγ` (§1378's `hζsphere`)

Combining §1130 (`‖ζ‖ ≤ M` on `ball (1+iγ) R''`) with §1397 (`sphere (σ+iγ) R ⊆ ball (1+iγ) R''` when
`σ-1+R < R''`) gives a ζ-bound on the sphere centred at the actual DLVP test point `σ+iγ` — exactly the `hζsphere`
hypothesis the two-radii factor bound §1378 needs.  EXIT0 / std-axioms.
-/

open Complex Metric

theorem exists_norm_riemannZeta_le_sphere {σ γ R R'' : ℝ} (hσ : 1 < σ) (hR''0 : 0 < R'')
    (hR''1 : R'' < 1) (hR''t : R'' < |γ|) (hincl : σ - 1 + R < R'') :
    ∃ Mζ : ℝ, ∀ s ∈ sphere ((σ : ℂ) + I * (γ : ℝ)) R, ‖riemannZeta s‖ ≤ Mζ := by
  obtain ⟨M, hM⟩ := exists_norm_riemannZeta_le_on_ball (t₀ := γ) (r := R'') hR''0 hR''1 hR''t
  refine ⟨M, fun s hs => hM s (sphere_subset_ball' ?_ hs)⟩
  have hdist : dist ((σ : ℂ) + I * (γ : ℝ)) (1 + I * (γ : ℝ)) = σ - 1 := by
    rw [Complex.dist_eq,
      show ((σ : ℂ) + I * (γ : ℝ)) - (1 + I * (γ : ℝ)) = (((σ - 1 : ℝ)) : ℂ) from by push_cast; ring,
      Complex.norm_real, Real.norm_eq_abs, abs_of_pos (by linarith)]
  rw [hdist]; exact hincl

#print axioms exists_norm_riemannZeta_le_sphere
