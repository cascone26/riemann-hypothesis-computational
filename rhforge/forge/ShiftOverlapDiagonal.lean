import ShiftOverlapIntegral
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic

/-!
# §1731 — Divisor-basis diagonal self-intersection (quantized shift-overlap)

Specializes §1730 to the actual interval-eigenbasis `vₖ(u) = sin(ωₖ u)` with the Dirichlet quantization
`ωₖ · L = k·π` (the boundary condition `vₖ(±a)=0`, `L = 2a`), on the geometrically-forced window `T = L−s`
(a divisor-bump at Arakelov distance `s = log pᵏ` overlaps its shifted copy only over `[0, L−s]`).  Under
quantization the boundary term collapses via `sin(2kπ − ωs) = −sin(ωs)`, giving the exact diagonal
self-intersection number:

`∫₀^{L−s} sin(ω(u+s))·sin(ω u) du = (L−s)/2 · cos(ω s) + sin(ω s)/(2ω)`.

Its symmetrised double `S_kk(s) = (L−s)cos(ωs) + sin(ωs)/ω` is exactly the `S_kk` used in the verified exact
Weil form (`suzuki_exact_Q`, §242).  This is the self-intersection entry of the finite Arakelov intersection
matrix — the diagonal the `ℓ²` toy model gets only in the `s→0` limit.  EXIT0 / std-axioms.
-/

open Real intervalIntegral

/-- **Quantized diagonal shift-overlap.**  For `ω ≠ 0` and `ω·L = k·π` (`k : ℤ`),
`∫₀^{L−s} sin(ω(u+s))·sin(ω u) du = (L−s)/2·cos(ω s) + sin(ω s)/(2ω)`. -/
theorem shift_overlap_diagonal (ω s L : ℝ) (k : ℤ) (hω : ω ≠ 0) (hL : ω * L = k * π) :
    (∫ u in (0:ℝ)..(L - s), Real.sin (ω * (u + s)) * Real.sin (ω * u))
      = (L - s) / 2 * Real.cos (ω * s) + Real.sin (ω * s) / (2 * ω) := by
  rw [shift_overlap_integral ω s (L - s) hω]
  have hkey : Real.sin (ω * (2 * (L - s) + s)) = - Real.sin (ω * s) := by
    have harg : ω * (2 * (L - s) + s) = -(ω * s) + (k : ℝ) * (2 * π) := by
      have hexp : ω * (2 * (L - s) + s) = 2 * (ω * L) - ω * s := by ring
      rw [hexp, hL]; push_cast; ring
    rw [harg, Real.sin_add_int_mul_two_pi, Real.sin_neg]
  rw [hkey]
  field_simp
  ring
