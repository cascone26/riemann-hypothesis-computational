import Mathlib.Analysis.SpecialFunctions.Trigonometric.Deriv
import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus

/-!
# §1732 — The off-diagonal shift-overlap (cross-prime intersection) integral in closed form

Companion to §1730: the **off-diagonal** entries of the finite Arakelov intersection matrix — the overlap of
two DISTINCT interval eigenmodes `sin(a·) , sin(b·)` shifted by the Arakelov distance `s`.  For `a ≠ ±b`,

`∫₀ᵀ sin(a(u+s))·sin(b u) du = (sin((a−b)T+a s)−sin(a s))/(2(a−b)) − (sin((a+b)T+a s)−sin(a s))/(2(a+b))`.

Together with §1730/§1731 this gives the ENTIRE finite intersection kernel `S_kl` of the exact Weil form in
closed form.  Under the Dirichlet quantization these off-diagonal entries carry the parity decoupling
(`k+l` odd ⟹ `S_kl=0`, verified numerically to 9e-16) that block-diagonalises the intersection form — the
next corollary.  Proof: FTC-2 against the antiderivative
`G(u) = sin((a−b)u+a s)/(2(a−b)) − sin((a+b)u+a s)/(2(a+b))`, product-to-sum `sin A sin B = (cos(A−B)−cos(A+B))/2`
discharging `G' =` integrand.  EXIT0 / std-axioms.
-/

open Real intervalIntegral

/-- **Off-diagonal shift-overlap, closed form.**  For `a − b ≠ 0` and `a + b ≠ 0`. -/
theorem shift_overlap_offdiag (a b s T : ℝ) (hd : a - b ≠ 0) (hs : a + b ≠ 0) :
    (∫ u in (0:ℝ)..T, Real.sin (a * (u + s)) * Real.sin (b * u))
      = (Real.sin ((a - b) * T + a * s) - Real.sin (a * s)) / (2 * (a - b))
        - (Real.sin ((a + b) * T + a * s) - Real.sin (a * s)) / (2 * (a + b)) := by
  have hderiv : ∀ u : ℝ, HasDerivAt
      (fun u : ℝ => Real.sin ((a - b) * u + a * s) / (2 * (a - b))
                      - Real.sin ((a + b) * u + a * s) / (2 * (a + b)))
      (Real.sin (a * (u + s)) * Real.sin (b * u)) u := by
    intro u
    have hg1 : HasDerivAt (fun u : ℝ => (a - b) * u + a * s) (a - b) u := by
      simpa using ((hasDerivAt_id u).const_mul (a - b)).add_const (a * s)
    have hg2 : HasDerivAt (fun u : ℝ => (a + b) * u + a * s) (a + b) u := by
      simpa using ((hasDerivAt_id u).const_mul (a + b)).add_const (a * s)
    have ht1 : HasDerivAt (fun u : ℝ => Real.sin ((a - b) * u + a * s) / (2 * (a - b)))
        (Real.cos ((a - b) * u + a * s) * (a - b) / (2 * (a - b))) u :=
      (hg1.sin).div_const (2 * (a - b))
    have ht2 : HasDerivAt (fun u : ℝ => Real.sin ((a + b) * u + a * s) / (2 * (a + b)))
        (Real.cos ((a + b) * u + a * s) * (a + b) / (2 * (a + b))) u :=
      (hg2.sin).div_const (2 * (a + b))
    have hGd := ht1.sub ht2
    -- product-to-sum value equality
    have hsub : Real.cos ((a - b) * u + a * s) =
        Real.cos (a * (u + s)) * Real.cos (b * u) + Real.sin (a * (u + s)) * Real.sin (b * u) := by
      have h := Real.cos_sub (a * (u + s)) (b * u)
      rw [show a * (u + s) - b * u = (a - b) * u + a * s by ring] at h
      exact h
    have hadd : Real.cos ((a + b) * u + a * s) =
        Real.cos (a * (u + s)) * Real.cos (b * u) - Real.sin (a * (u + s)) * Real.sin (b * u) := by
      have h := Real.cos_add (a * (u + s)) (b * u)
      rw [show a * (u + s) + b * u = (a + b) * u + a * s by ring] at h
      exact h
    have hval : Real.sin (a * (u + s)) * Real.sin (b * u)
        = Real.cos ((a - b) * u + a * s) * (a - b) / (2 * (a - b))
          - Real.cos ((a + b) * u + a * s) * (a + b) / (2 * (a + b)) := by
      have e1 : Real.cos ((a - b) * u + a * s) * (a - b) / (2 * (a - b))
          = Real.cos ((a - b) * u + a * s) / 2 := by field_simp
      have e2 : Real.cos ((a + b) * u + a * s) * (a + b) / (2 * (a + b))
          = Real.cos ((a + b) * u + a * s) / 2 := by field_simp
      rw [e1, e2, hsub, hadd]; ring
    rw [hval]; exact hGd
  have hcont : Continuous (fun u : ℝ => Real.sin (a * (u + s)) * Real.sin (b * u)) := by
    fun_prop
  rw [integral_eq_sub_of_hasDerivAt (fun u _ => hderiv u) (hcont.intervalIntegrable 0 T)]
  simp only [mul_zero, zero_add]
  ring
