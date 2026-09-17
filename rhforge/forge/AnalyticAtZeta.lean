import Mathlib.NumberTheory.LSeries.RiemannZeta
import Mathlib.Analysis.Complex.CauchyIntegral

/-!
# §1345 — ζ is analytic away from its pole: `z ≠ 1 ⟹ AnalyticAt ℂ ζ z`

ζ is differentiable on the open set `{1}ᶜ` (`differentiableAt_riemannZeta`), so by `DifferentiableOn.analyticAt`
it is analytic at every `z ≠ 1`.  Combined with §1344 this gives `ContinuousAt (logDeriv ζ) z` at the DLVP test
point `σ+iγ` (`Re > 1`, `≠ 1`, `ζ ≠ 0`).  EXIT0 / std-axioms.
-/

open Complex Topology

theorem analyticAt_riemannZeta {z : ℂ} (hz : z ≠ 1) : AnalyticAt ℂ riemannZeta z := by
  refine DifferentiableOn.analyticAt (s := {(1 : ℂ)}ᶜ)
    (fun w hw => (differentiableAt_riemannZeta hw).differentiableWithinAt) ?_
  exact isOpen_compl_singleton.mem_nhds hz

#print axioms analyticAt_riemannZeta
