import Mathlib.Analysis.SpecialFunctions.Complex.Log
import Mathlib.Analysis.Calculus.InverseFunctionTheorem.Deriv
import Mathlib.Analysis.SpecialFunctions.Complex.Analytic

/-!
# §1323 — A continuous log-branch of a holomorphic function is holomorphic

If `φ` is continuous at `z₀` and `exp∘φ = ψ` near `z₀` with `ψ` differentiable at `z₀`, then `φ` is
differentiable at `z₀` with `φ'(z₀) = ψ'(z₀)/ψ(z₀) = logDeriv ψ z₀`.  Proof via the inverse function theorem for
`exp` (nonzero derivative): let `L` be the local inverse of `exp` near `w₀ = φ z₀`; then `φ = L∘ψ` near `z₀`
(continuity + `L∘exp = id`), and the chain rule gives the derivative.  This is the continuous→holomorphic UPGRADE
of the branch §1322 — the last piece: with §1320/§1321/§1318 it yields the Borel–Carathéodory `Re(logDeriv(ζ/P)) ≥ -B`.
EXIT0 / std-axioms.
-/

open Complex Filter Topology

theorem hasDerivAt_of_exp_comp {φ ψ : ℂ → ℂ} {z₀ dψ : ℂ}
    (hφ : ContinuousAt φ z₀) (hψ : HasDerivAt ψ dψ z₀)
    (hcomp : ∀ᶠ z in 𝓝 z₀, Complex.exp (φ z) = ψ z) :
    HasDerivAt φ (dψ / ψ z₀) z₀ := by
  have hexp0 : Complex.exp (φ z₀) = ψ z₀ := hcomp.self_of_nhds
  have hf' : Complex.exp (φ z₀) ≠ 0 := Complex.exp_ne_zero _
  have hsd : HasStrictDerivAt Complex.exp (Complex.exp (φ z₀)) (φ z₀) := Complex.hasStrictDerivAt_exp _
  set L := hsd.localInverse Complex.exp (Complex.exp (φ z₀)) (φ z₀) hf' with hLdef
  have hLd : HasStrictDerivAt L (Complex.exp (φ z₀))⁻¹ (Complex.exp (φ z₀)) := hsd.to_localInverse hf'
  have hLexp : ∀ᶠ x in 𝓝 (φ z₀), L (Complex.exp x) = x := hsd.eventually_left_inverse hf'
  have hφLψ : φ =ᶠ[𝓝 z₀] (fun z => L (ψ z)) := by
    have h1 : ∀ᶠ z in 𝓝 z₀, L (Complex.exp (φ z)) = φ z := hφ.eventually hLexp
    filter_upwards [h1, hcomp] with z hz hc
    rw [← hc, hz]
  have hLd' : HasDerivAt L (Complex.exp (φ z₀))⁻¹ (ψ z₀) := hexp0 ▸ hLd.hasDerivAt
  have hcomp2 : HasDerivAt (fun z => L (ψ z)) ((Complex.exp (φ z₀))⁻¹ * dψ) z₀ :=
    hLd'.comp z₀ hψ
  have hval : dψ / ψ z₀ = (Complex.exp (φ z₀))⁻¹ * dψ := by
    rw [hexp0, div_eq_mul_inv, mul_comm]
  rw [hval]
  exact hcomp2.congr_of_eventuallyEq hφLψ

#print axioms hasDerivAt_of_exp_comp
