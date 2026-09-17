import ZetaStripBound

set_option maxHeartbeats 1000000

open Complex Metric

/-!
# §1684 — the CONCRETE sphere bound at the Landau center: `‖ζ‖ ≤ 5T` on `sphere(1+ε+iT, 1/4)`

Instantiates the 8/12 strip growth bound `‖ζ(s)‖ ≤ ‖s‖/‖s-1‖ + ‖s‖/Re s` on the Landau ball's
boundary: for the center `c = 1+ε+iT` (`0 < ε ≤ ½`, `T ≥ 2`), every sphere point has
`Re ≥ ¾`, `‖s‖ ≤ 2T`, `‖s-1‖ ≥ T/2`, giving `‖ζ‖ ≤ 4 + 8T/3 ≤ 5T`. The `M₀`-input of §1683:
`log((M₀+1)/‖ζ c‖) = O(log T)`.
-/

theorem zeta_sphere_bound_center {T ε : ℝ} (hT : 2 ≤ T) (hε : 0 < ε) (hε' : ε ≤ 1 / 2) :
    ∀ s ∈ sphere (((1 + ε : ℝ) : ℂ) + (T : ℝ) * Complex.I) (1 / 4),
      ‖riemannZeta s‖ ≤ 5 * T := by
  intro s hs
  set c : ℂ := ((1 + ε : ℝ) : ℂ) + (T : ℝ) * Complex.I with hcdef
  have hsc : ‖s - c‖ = 1 / 4 := by
    rw [← dist_eq_norm]
    exact mem_sphere.mp hs
  have hcre : c.re = 1 + ε := by rw [hcdef]; simp
  have hcim : c.im = T := by rw [hcdef]; simp
  -- real part of s
  have hre_diff : |(s - c).re| ≤ 1 / 4 := by
    calc |(s - c).re| ≤ ‖s - c‖ := abs_re_le_norm _
      _ = 1 / 4 := hsc
  have hsre : 3 / 4 ≤ s.re := by
    have h1 : (s - c).re = s.re - c.re := by simp
    have := abs_le.mp hre_diff
    rw [h1, hcre] at this
    linarith [this.1]
  -- imaginary part of s
  have him_diff : |(s - c).im| ≤ 1 / 4 := by
    calc |(s - c).im| ≤ ‖s - c‖ := abs_im_le_norm _
      _ = 1 / 4 := hsc
  have hsim : T - 1 / 4 ≤ s.im := by
    have h1 : (s - c).im = s.im - c.im := by simp
    have := abs_le.mp him_diff
    rw [h1, hcim] at this
    linarith [this.1]
  have hs1 : s ≠ 1 := by
    intro h
    rw [h] at hsim
    simp at hsim
    linarith
  -- norm bounds
  have hcnorm : ‖c‖ ≤ 1 + ε + T := by
    rw [hcdef]
    calc ‖((1 + ε : ℝ) : ℂ) + (T : ℝ) * Complex.I‖
        ≤ ‖((1 + ε : ℝ) : ℂ)‖ + ‖(T : ℝ) * Complex.I‖ := norm_add_le _ _
      _ = |1 + ε| + |T| := by
          rw [Complex.norm_real, Real.norm_eq_abs, norm_mul, Complex.norm_real,
            Real.norm_eq_abs, Complex.norm_I, mul_one]
      _ = 1 + ε + T := by
          rw [abs_of_pos (by linarith), abs_of_pos (by linarith)]
  have hsnorm : ‖s‖ ≤ 2 * T := by
    calc ‖s‖ = ‖(s - c) + c‖ := by ring_nf
      _ ≤ ‖s - c‖ + ‖c‖ := norm_add_le _ _
      _ ≤ 1 / 4 + (1 + ε + T) := by rw [hsc]; linarith
      _ ≤ 2 * T := by linarith
  have hs1norm : T / 2 ≤ ‖s - 1‖ := by
    calc T / 2 ≤ T - 1 / 4 := by linarith
      _ ≤ s.im := hsim
      _ ≤ |(s - 1).im| := by
          rw [show (s - 1 : ℂ).im = s.im by simp]
          exact le_abs_self _
      _ ≤ ‖s - 1‖ := abs_im_le_norm _
  -- assemble via the strip bound
  have hstrip := norm_riemannZeta_le_strip (lt_of_lt_of_le (by norm_num) hsre) hs1
  have hd1 : ‖s‖ / ‖s - 1‖ ≤ 2 * T / (T / 2) := by
    gcongr <;> linarith
  have hd2 : ‖s‖ / s.re ≤ 2 * T / (3 / 4) := by
    gcongr <;> linarith
  have he1 : 2 * T / (T / 2) = 4 := by
    field_simp
    ring
  have he2 : 2 * T / (3 / 4) = 8 * T / 3 := by
    ring
  calc ‖riemannZeta s‖ ≤ ‖s‖ / ‖s - 1‖ + ‖s‖ / s.re := hstrip
    _ ≤ 2 * T / (T / 2) + 2 * T / (3 / 4) := by linarith
    _ = 4 + 8 * T / 3 := by rw [he1, he2]
    _ ≤ 5 * T := by linarith

#print axioms zeta_sphere_bound_center
