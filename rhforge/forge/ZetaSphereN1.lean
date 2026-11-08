import ZetaSphereUniform

/-!
# §1447 — simplified sphere bound at `N = 1`

Specializing §1440 to `N = 1` collapses every `1^exponent = 1`, killing the `rpow` terms and leaving the clean
explicit bound

  `‖ζ z‖ ≤ 2 + 1/(‖c-1‖ - R) + (‖c‖ + R)/(c.re - R)`   for `z ∈ sphere c R`, `0 < c.re - R`, `R < ‖c-1‖`.

Still `≍ |γ|` on the DLVP sphere (`‖c‖ ≍ |γ|`, `c.re - R ≍ 1`), but with no exponentials — the version the final
constant-choice assembly threads.  EXIT0 / std-axioms.
-/

open Metric

theorem norm_riemannZeta_le_sphere_N1 {c : ℂ} {R : ℝ} (hR : 0 < R) (hcR : 0 < c.re - R)
    (hc1 : R < ‖c - 1‖) {z : ℂ} (hz : z ∈ sphere c R) :
    ‖riemannZeta z‖ ≤ 2 + 1 / (‖c - 1‖ - R) + (‖c‖ + R) / (c.re - R) := by
  have h := norm_riemannZeta_le_sphere_uniform (N := 1) le_rfl hR hcR hc1 hz
  simp only [Nat.cast_one, Real.one_rpow] at h
  calc ‖riemannZeta z‖
      ≤ (1 : ℝ) + 1 + 1 / (‖c - 1‖ - R) + (‖c‖ + R) * (1 / (c.re - R)) := h
    _ = 2 + 1 / (‖c - 1‖ - R) + (‖c‖ + R) / (c.re - R) := by ring

#print axioms norm_riemannZeta_le_sphere_N1
