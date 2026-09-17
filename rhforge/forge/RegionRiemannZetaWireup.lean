import RegionZpow
import DLVPPositivityZeta

/-!
# §1406 — ζ region wire-up: `hdlvp` auto-discharged, abstract in the factorization data

Applies the region master §1336 at the DLVP test point `σ_K+iγ` (`σ_K = 1+1/(2K)`, `K = 3A+4B+C`), discharging the
`3-4-1` positivity `hdlvp` internally via §1358 (`dlvp_positivity_riemannZeta`).  Abstract in the zero-set `S`,
multiplicities `m`, factor value `g`, and the decomposition/BC/pole/growth bounds — so the final assembly need only
supply those from §1400 (`hdecomp`), §1402 (`hgB`), §1405 (`hpole`), §1404 (`hgrowth`) plus the geometric
multiplicity facts.  Conclusion: `β ≤ 1 - 1/(14K)`.  This is the quantitative de la Vallée Poussin zero-free region
for ζ, modulo the constant-choice fixed-point.  EXIT0 / std-axioms.
-/

open Complex

theorem region_riemannZeta_wireup {γ β A B C : ℝ} {g : ℂ} {S : Finset ℂ} {m : ℂ → ℤ}
    (hK : 0 < 3 * A + 4 * B + C) (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hβ : β < 1 + 1 / (2 * (3 * A + 4 * B + C)))
    (hρ₀ : ((β : ℂ) + I * (γ : ℝ)) ∈ S)
    (hleft : ∀ ρ ∈ S, ρ.re < 1 + 1 / (2 * (3 * A + 4 * B + C)))
    (hm : ∀ ρ ∈ S, 0 ≤ m ρ) (hm₀ : 1 ≤ m ((β : ℂ) + I * (γ : ℝ)))
    (hgB : -B ≤ g.re)
    (hdecomp : logDeriv riemannZeta (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ))
      = (∑ ρ ∈ S, (m ρ : ℂ) / (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ) - ρ)) + g)
    (hpole : (-logDeriv riemannZeta (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ))).re
      ≤ 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1) + A)
    (hgrowth : (-logDeriv riemannZeta
        (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ))).re ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hσK : (1 : ℝ) < 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
    have : (0 : ℝ) < 1 / (2 * (3 * A + 4 * B + C)) := by positivity
    linarith
  have hdlvp := dlvp_positivity_riemannZeta (σ := 1 + 1 / (2 * (3 * A + 4 * B + C))) (γ := γ) hσK
  exact region_of_zpow_decomp_and_bounds hK hA hB hC hβ hρ₀ hleft hm hm₀ hgB hdecomp hdlvp hpole hgrowth

#print axioms region_riemannZeta_wireup
