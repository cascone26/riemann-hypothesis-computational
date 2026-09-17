import LFunctionSphereUniform

/-!
# §1535 — scalar sphere bound for `LFunction χ` (`N = 1` collapse)

L-analogue of ζ §1447 (`norm_riemannZeta_le_sphere_N1`).  Specializes the sphere-uniform bound §1534 at `N = 1`,
which collapses the `(N:ℝ)^{…}` factor to `1` (`Real.one_rpow`) and turns `((N+1:ℕ):ℝ)` into `2`.  Gives a clean
purely-scalar bound on `‖LFunction χ‖` over the sphere with NO `rpow` on a variable base except the fixed `2^{…}`:

  `‖LFunction χ z‖ ≤ 1 + (M · 2^{-(c.re-R)} + (‖c‖+R) · M / (c.re-R))`   (`z ∈ sphere c R`).

`M` is the character modulus.  With `c = σK+iγ`, `R = 1/4`: `Ms ≍ M·|γ|`, so `log Ms ≍ log|γ|` — the input the
sharp fixed-point consumes.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric

theorem norm_LFunction_le_sphere_N1 {M : ℕ} [NeZero M] {χ : DirichletCharacter ℂ M} (hχ : χ ≠ 1)
    {c : ℂ} {R : ℝ} (hR : 0 < R) (hcR : 0 < c.re - R) {z : ℂ} (hz : z ∈ sphere c R) :
    ‖LFunction χ z‖ ≤ 1
      + ((M : ℝ) * (2 : ℝ) ^ (-(c.re - R)) + (‖c‖ + R) * M / (c.re - R)) := by
  have h := norm_LFunction_le_sphere_uniform hχ (le_refl 1) hR hcR hz
  simp only [Nat.cast_one, Real.one_rpow, mul_one] at h
  calc ‖LFunction χ z‖
      ≤ 1 + ((2 : ℝ) ^ (-(c.re - R)) * M + (‖c‖ + R) / (c.re - R) * M) := h
    _ = 1 + ((M : ℝ) * (2 : ℝ) ^ (-(c.re - R)) + (‖c‖ + R) * M / (c.re - R)) := by ring

#print axioms norm_LFunction_le_sphere_N1
