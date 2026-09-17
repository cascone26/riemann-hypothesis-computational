import WindingCoreFTC

set_option maxHeartbeats 1000000

open Complex intervalIntegral

/-!
# §1655 — WINDING EDGES: closed forms for the four rectangle-edge integrals of `(z - p)⁻¹`

Both horizontal edges (`z = x + c·I`) and vertical edges (`z = a + y·I`) of a rectangle avoiding
`p` reduce to the §1654 `windingPrim` closed form. Rather than substitution lemmas (whose cast
matching is fragile), each edge gets a direct FTC with the SHIFTED antiderivative
`x ↦ windingPrim k (x - p.re)` (resp. `y ↦ -I·windingPrim k (y - p.im)`) — the vertical case
absorbing the factor from `h + v·I = (v + (-h)·I)·I`.
-/

/-- §1655a — horizontal edge: `∫_{x₁}^{x₂} ((x + c·I) - p)⁻¹ dx` in closed `windingPrim` form,
for any horizontal line `Im = c` avoiding `p`. -/
theorem integral_horizontal_winding (p : ℂ) {c : ℝ} (hc : c ≠ p.im) (x₁ x₂ : ℝ) :
    ∫ x in x₁..x₂, ((x : ℂ) + c * Complex.I - p)⁻¹
      = windingPrim (c - p.im) (x₂ - p.re) - windingPrim (c - p.im) (x₁ - p.re) := by
  have hk : c - p.im ≠ 0 := sub_ne_zero.mpr hc
  have hne : ∀ x : ℝ, (x : ℂ) + c * Complex.I - p ≠ 0 := by
    intro x h
    have him := congrArg Complex.im h
    simp at him
    exact hc (by linarith)
  have hderiv : ∀ x : ℝ, HasDerivAt (fun x : ℝ => windingPrim (c - p.im) (x - p.re))
      (((x : ℂ) + c * Complex.I - p)⁻¹) x := by
    intro x
    have hshift : HasDerivAt (fun x : ℝ => x - p.re) 1 x := (hasDerivAt_id x).sub_const p.re
    have h1 := (hasDerivAt_windingPrim hk (x - p.re)).scomp x hshift
    have hval : ((x : ℂ) + c * Complex.I - p)⁻¹
        = ((↑(x - p.re) : ℂ) + ↑(c - p.im) * Complex.I)⁻¹ := by
      congr 1
      push_cast
      linear_combination Complex.re_add_im p
    rw [hval]
    simpa [Function.comp_def] using h1
  have hcont : Continuous fun x : ℝ => ((x : ℂ) + c * Complex.I - p)⁻¹ :=
    Continuous.inv₀ ((Complex.continuous_ofReal.add continuous_const).sub continuous_const) hne
  exact integral_eq_sub_of_hasDerivAt (fun x _ => hderiv x) (hcont.intervalIntegrable x₁ x₂)

/-- §1655b — vertical edge: `∫_{y₁}^{y₂} ((a + y·I) - p)⁻¹ dy` in closed `windingPrim` form,
for any vertical line `Re = a` avoiding `p`. The `-I` factor comes from
`h + v·I = (v + (-h)·I)·I`. -/
theorem integral_vertical_winding (p : ℂ) {a : ℝ} (ha : a ≠ p.re) (y₁ y₂ : ℝ) :
    ∫ y in y₁..y₂, ((a : ℂ) + y * Complex.I - p)⁻¹
      = -Complex.I * windingPrim (p.re - a) (y₂ - p.im)
        - -Complex.I * windingPrim (p.re - a) (y₁ - p.im) := by
  have hk : p.re - a ≠ 0 := sub_ne_zero.mpr (Ne.symm ha)
  have hne : ∀ y : ℝ, (a : ℂ) + y * Complex.I - p ≠ 0 := by
    intro y h
    have hre := congrArg Complex.re h
    simp at hre
    exact ha (by linarith)
  have hbase : ∀ y : ℝ, (a : ℂ) + y * Complex.I - p
      = ((↑(y - p.im) : ℂ) + ↑(p.re - a) * Complex.I) * Complex.I := by
    intro y
    push_cast
    linear_combination Complex.re_add_im p + ((a : ℂ) - ↑p.re) * Complex.I_sq
  have hderiv : ∀ y : ℝ, HasDerivAt (fun y : ℝ => -Complex.I * windingPrim (p.re - a) (y - p.im))
      (((a : ℂ) + y * Complex.I - p)⁻¹) y := by
    intro y
    have hshift : HasDerivAt (fun y : ℝ => y - p.im) 1 y := (hasDerivAt_id y).sub_const p.im
    have h1 := ((hasDerivAt_windingPrim hk (y - p.im)).scomp y hshift).const_mul (-Complex.I)
    have hval : ((a : ℂ) + y * Complex.I - p)⁻¹
        = -Complex.I * ((↑(y - p.im) : ℂ) + ↑(p.re - a) * Complex.I)⁻¹ := by
      rw [hbase y, mul_inv_rev, Complex.inv_I]
    rw [hval]
    simpa [Function.comp_def] using h1
  have hcont : Continuous fun y : ℝ => ((a : ℂ) + y * Complex.I - p)⁻¹ :=
    Continuous.inv₀
      ((continuous_const.add ((Complex.continuous_ofReal.mul continuous_const))).sub
        continuous_const) hne
  exact integral_eq_sub_of_hasDerivAt (fun y _ => hderiv y) (hcont.intervalIntegrable y₁ y₂)

#print axioms integral_horizontal_winding
#print axioms integral_vertical_winding
