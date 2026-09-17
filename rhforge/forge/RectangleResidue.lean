import RectangleWinding

set_option maxHeartbeats 1000000

open Complex intervalIntegral

/-!
# §1658 — RECTANGLE RESIDUE (simple pole): `∮_∂R c·(z - p)⁻¹ dz = 2πi·c` for `p` inside `R`

The residue-scaled form of §1657: a simple pole with residue `c` deposits `2πi·c` on any rectangle
strictly containing it. Applied with `c = x²/2`, `p = 1`, this is EXACTLY the deposit the Perron
pole term `(x²/2)/(s-1)` makes when the Perron contour shifts left past `s = 1` — the main-term
mechanism of PNT-with-error. Proof: pull the constant out of all four §1657 edge integrals.
-/

theorem integral_boundary_rect_residue (c p z w : ℂ)
    (h1 : z.re < p.re) (h2 : p.re < w.re) (h3 : z.im < p.im) (h4 : p.im < w.im) :
    (∫ x in z.re..w.re, c * ((x : ℂ) + z.im * Complex.I - p)⁻¹) -
      (∫ x in z.re..w.re, c * ((x : ℂ) + w.im * Complex.I - p)⁻¹) +
      Complex.I • (∫ y in z.im..w.im, c * ((w.re : ℂ) + y * Complex.I - p)⁻¹) -
      Complex.I • (∫ y in z.im..w.im, c * ((z.re : ℂ) + y * Complex.I - p)⁻¹)
      = 2 * (Real.pi : ℂ) * Complex.I * c := by
  have hbase := integral_boundary_rect_sub_inv p z w h1 h2 h3 h4
  simp only [integral_const_mul, smul_eq_mul] at *
  linear_combination c * hbase

#print axioms integral_boundary_rect_residue
