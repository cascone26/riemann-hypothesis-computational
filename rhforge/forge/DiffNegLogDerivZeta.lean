import Mathlib.NumberTheory.LSeries.RiemannZeta
import Mathlib.Analysis.Complex.CauchyIntegral
import Mathlib.Analysis.Calculus.FDeriv.Analytic
import Mathlib.Analysis.Calculus.LogDeriv

set_option maxHeartbeats 1000000

/-!
# §1606 — holomorphy of `-ζ'/ζ` on the zero-free-and-not-`1` set (contour-shift differentiability input)

`-ζ'/ζ(s) = -(deriv ζ / ζ)(s)` is complex-differentiable at every `s ≠ 1` with `ζ s ≠ 0`: `ζ` is analytic on
`{1}ᶜ` (`analyticOn_riemannZeta`), so both `ζ` and `ζ'` are differentiable there, and `ζ s ≠ 0` lets the quotient
differentiate.  This is the pointwise holomorphy needed to feed the Perron integrand
`s ↦ -ζ'/ζ(s)·x^{s+1}/(s(s+1))` into Mathlib's rectangle-contour vanishing lemma
(`integral_boundary_rect_eq_zero_of_differentiable_on_off_countable`) after subtracting the `s=1` principal part
(residue §1603 = `x²/2`).  EXIT0 / std-axioms.
-/

open Complex

theorem differentiableAt_neg_logDeriv_riemannZeta {s : ℂ} (hs1 : s ≠ 1)
    (hζ : riemannZeta s ≠ 0) :
    DifferentiableAt ℂ (fun z => -logDeriv riemannZeta z) s := by
  have hmem : s ∈ ({1}ᶜ : Set ℂ) := hs1
  have hderiv_diff : DifferentiableAt ℂ (deriv riemannZeta) s :=
    (analyticOn_riemannZeta.deriv s hmem).differentiableAt
  have hζ_diff : DifferentiableAt ℂ riemannZeta s := differentiableAt_riemannZeta hs1
  have hfun : (fun z => -logDeriv riemannZeta z)
      = (fun z => -(deriv riemannZeta z / riemannZeta z)) := by
    funext z; rw [logDeriv_apply]
  rw [hfun]
  exact (hderiv_diff.div hζ_diff hζ).neg

#print axioms differentiableAt_neg_logDeriv_riemannZeta
