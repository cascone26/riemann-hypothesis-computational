import Mathlib.Analysis.SpecialFunctions.Trigonometric.Deriv
import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus

/-!
# §1730 — The shift-overlap (finite intersection) integral in closed form

The finite-place content of the Weil explicit-formula pairing (verified 2026-08-24 on the exact Suzuki form
`suzuki_exact_Q.build_Q`: the prime piece is `−Σ Λ(n)/√n · S_kl(log n)`, coefficient exactly `Λ(pᵏ)/p^{k/2}`)
reduces to the **shift-overlap kernel** `S_kl(s) = ∫ vₖ(y+s) vₗ(y) dy` of interval basis functions — the local
intersection multiplicity of two divisor-bumps at Arakelov "distance" `s = log pᵏ`.  Here we bank its exact
closed form for the diagonal generator `vₖ(u) = sin(ω u)`:

`∫₀ᵀ sin(ω(u+s)) · sin(ω u) du = (T/2) cos(ω s) − (sin(ω(2T+s)) − sin(ω s))/(4ω)`  for `ω ≠ 0`.

This is the exact self-intersection number that upgrades the Lean Arakelov pairing from its `ℓ²` (`δ_kl`) toy
to the true translation-overlap.  Proof: FTC-2 against the explicit antiderivative
`F(u) = (cos(ω s)/2)·u − sin(ω(2u+s))/(4ω)`, with the product-to-sum identity
`sin A · sin B = (cos(A−B) − cos(A+B))/2` discharging `F' = ` integrand.  EXIT0 / std-axioms.
-/

open Real intervalIntegral

/-- **Shift-overlap integral, closed form.**  For `ω ≠ 0`,
`∫₀ᵀ sin(ω(u+s))·sin(ω u) du = (T/2)cos(ω s) − (sin(ω(2T+s)) − sin(ω s))/(4ω)`. -/
theorem shift_overlap_integral (ω s T : ℝ) (hω : ω ≠ 0) :
    (∫ u in (0:ℝ)..T, Real.sin (ω * (u + s)) * Real.sin (ω * u))
      = T / 2 * Real.cos (ω * s)
        - (Real.sin (ω * (2 * T + s)) - Real.sin (ω * s)) / (4 * ω) := by
  have hderiv : ∀ u : ℝ, HasDerivAt
      (fun u : ℝ => Real.cos (ω * s) / 2 * u - Real.sin (ω * (2 * u + s)) / (4 * ω))
      (Real.sin (ω * (u + s)) * Real.sin (ω * u)) u := by
    intro u
    -- inner map g u = ω*(2*u+s), g' = ω*2
    have hg : HasDerivAt (fun u : ℝ => ω * (2 * u + s)) (ω * 2) u := by
      have h1 : HasDerivAt (fun u : ℝ => 2 * u + s) (2 * 1 : ℝ) u :=
        (((hasDerivAt_id u).const_mul (2 : ℝ))).add_const s
      simpa using h1.const_mul ω
    have ht1 : HasDerivAt (fun u : ℝ => Real.cos (ω * s) / 2 * u) (Real.cos (ω * s) / 2) u := by
      simpa using (hasDerivAt_id u).const_mul (Real.cos (ω * s) / 2)
    have ht2 : HasDerivAt (fun u : ℝ => Real.sin (ω * (2 * u + s)) / (4 * ω))
        (Real.cos (ω * (2 * u + s)) * (ω * 2) / (4 * ω)) u :=
      (hg.sin).div_const (4 * ω)
    have hFd := ht1.sub ht2
    -- product-to-sum value equality
    have hsub : Real.cos (ω * s) =
        Real.cos (ω * (u + s)) * Real.cos (ω * u) + Real.sin (ω * (u + s)) * Real.sin (ω * u) := by
      have h := Real.cos_sub (ω * (u + s)) (ω * u)
      rw [show ω * (u + s) - ω * u = ω * s by ring] at h
      exact h
    have hadd : Real.cos (ω * (2 * u + s)) =
        Real.cos (ω * (u + s)) * Real.cos (ω * u) - Real.sin (ω * (u + s)) * Real.sin (ω * u) := by
      have h := Real.cos_add (ω * (u + s)) (ω * u)
      rw [show ω * (u + s) + ω * u = ω * (2 * u + s) by ring] at h
      exact h
    have hval : Real.sin (ω * (u + s)) * Real.sin (ω * u)
        = Real.cos (ω * s) / 2 - Real.cos (ω * (2 * u + s)) * (ω * 2) / (4 * ω) := by
      have h2 : Real.cos (ω * (2 * u + s)) * (ω * 2) / (4 * ω)
          = Real.cos (ω * (2 * u + s)) / 2 := by
        field_simp; ring
      rw [h2, hsub, hadd]; ring
    rw [hval]; exact hFd
  have hcont : Continuous (fun u : ℝ => Real.sin (ω * (u + s)) * Real.sin (ω * u)) := by
    fun_prop
  rw [integral_eq_sub_of_hasDerivAt (fun u _ => hderiv u) (hcont.intervalIntegrable 0 T)]
  simp only [mul_zero, zero_add]
  ring
