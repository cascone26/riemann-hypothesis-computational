import Mathlib.Data.Complex.Basic

/-!
# §1569 — `Re(1/z) ≤ 1/Re z` for `Re z > 0`

Opening brick of the SMALL-|t| / exceptional-zero arc.  `Re(1/z) = Re z / |z|² ≤ Re z / (Re z)² = 1/Re z`.
This is the complex pole-term bound the small-|t| region needs: the T₂ factor's pole contribution at
`s = σ+2iγ` obeys `Re(1/(s-1)) ≤ 1/(σ-1)` — same ceiling as the real-axis pole bound, uniformly in `γ`.
EXIT0 / std-axioms.
-/

open Complex

theorem re_inv_le_inv_re {z : ℂ} (hz : 0 < z.re) : (z⁻¹).re ≤ 1 / z.re := by
  rw [Complex.inv_re]
  have hsq : z.re * z.re ≤ Complex.normSq z := by
    rw [Complex.normSq_apply]
    exact le_add_of_nonneg_right (mul_self_nonneg z.im)
  have hnpos : 0 < Complex.normSq z := lt_of_lt_of_le (mul_pos hz hz) hsq
  rw [div_le_div_iff₀ hnpos hz, one_mul]
  exact hsq

#print axioms re_inv_le_inv_re
