import Psi1CMellinConvergent
import VerticalIntegrablePsi1C
import Mathlib.Analysis.MellinInversion

set_option maxHeartbeats 1000000

/-!
# §1601 ★★★★★ — THE SMOOTHED PERRON FORMULA: `Ψ(x) = mellinInv (-σ₀-1) (mellin Ψ) x` (`x > 0`, `σ₀ > 1`)

The SECOND PILLAR of the explicit formula.  All three hypotheses of Mathlib's `mellinInv_mellin_eq` are now proven:
`MellinConvergent Ψ (-σ₀-1)` (§1599), `VerticalIntegrable (mellin Ψ) (-σ₀-1)` (§1600), and `ContinuousAt Ψ x`
(free, §1590).  Unfolding `mellinInv` (and using §1598) gives the smoothed Perron representation

  `Ψ(x) = (1/2π)·∫_ℝ x^{-(σ_M+iy)}·(-ζ'/ζ(s_y))/(s_y(s_y+1)) dy`,

recovering the smoothed Chebyshev function as a vertical integral of `-ζ'/ζ`.  Contour-shifting this against the
zero-free region §1575 yields `Ψ(x) = x²/2 + O(error)`, and de-smoothing yields `ψ(x) = x + O(x·exp(-c√log x))`
— the prime number theorem with error term.  EXIT0 / std-axioms.
-/

open Complex Chebyshev MeasureTheory

theorem smoothed_perron {σ₀ x : ℝ} (hσ₀ : 1 < σ₀) (hx : 0 < x) :
    mellinInv (-σ₀ - 1) (mellin psi1C) x = psi1C x :=
  mellinInv_mellin_eq (-σ₀ - 1) psi1C hx (mellinConvergent_psi1C hσ₀)
    (verticalIntegrable_mellin_psi1C hσ₀) continuous_psi1C.continuousAt

#print axioms smoothed_perron
