import LogDerivEqCodiscretePt
import SmulFactorizedEval
import LogDerivZpowMulSplit

/-!
# §1342 — `logDeriv ζ` decomposition from the extract_zeros_poles codiscrete factorization

Threading converter #3 (§1341) + converter #2 (§1338) + the multiplicity split (§1332).  Given
`ζ =ᶠ[codiscreteWithin U] (∏ᶠ u,(·-u)^{d u})•g` (the extract_zeros_poles output, `d` the divisor, `g` nonvanishing
analytic) and, at a test point `z ∈ U` where both log-derivatives are continuous, `z` avoids the divisor support,
`g z ≠ 0`, `g` differentiable at `z`:

  `logDeriv ζ z = ∑_{ρ ∈ support d} d ρ/(z-ρ) + logDeriv g z`.

This is exactly the `hdecomp` §1335/§1336 consume — the ζ'/ζ decomposition obtained from Mathlib's divisor
factorization, multiplicity-general.  EXIT0 / std-axioms.
-/

open Filter Topology Function

theorem logDeriv_eq_sum_of_codiscrete_factor {ζ g : ℂ → ℂ} {d : ℂ → ℤ} {U : Set ℂ} {z : ℂ}
    (hU : IsOpen U) (hz : z ∈ U) (hdf : (Function.support d).Finite)
    (h : ζ =ᶠ[codiscreteWithin U] (∏ᶠ u, (· - u) ^ (d u)) • g)
    (hζc : ContinuousAt (logDeriv ζ) z)
    (hFc : ContinuousAt (logDeriv ((∏ᶠ u, (· - u) ^ (d u)) • g)) z)
    (hsub : ∀ ρ ∈ hdf.toFinset, z - ρ ≠ 0) (hgz : g z ≠ 0) (hgd : DifferentiableAt ℂ g z) :
    logDeriv ζ z = (∑ ρ ∈ hdf.toFinset, (d ρ : ℂ) / (z - ρ)) + logDeriv g z := by
  rw [logDeriv_eq_of_codiscreteWithin_of_continuousAt hU hz h hζc hFc]
  have hFeq : ((∏ᶠ u, (· - u) ^ (d u)) • g) = fun s => (∏ ρ ∈ hdf.toFinset, (s - ρ) ^ (d ρ)) * g s := by
    funext s; exact smul_factorized_rational_eval hdf s
  rw [hFeq]
  exact logDeriv_prod_zpow_mul_split hsub hgz hgd

#print axioms logDeriv_eq_sum_of_codiscrete_factor
