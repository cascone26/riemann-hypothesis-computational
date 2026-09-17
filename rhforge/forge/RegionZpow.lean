import HzeroZpow
import DLVPRegionFromBounds

/-!
# §1336 ★ Multiplicity master composition: region from the ζ decomposition (any multiplicity) + 3 bounds

The multiplicity version of §1326: chaining §1335 (`Re(-F) ≤ B - 1/(σ-β)` from the divisor-weighted decomposition)
into §1311 (four bounds ⟹ region).  At `σ = 1 + 1/(2K)`, `K = 3A+4B+C`, given the `3-4-1` positivity with
`T₁ = Re(-F)` (`F = logDeriv ζ(σ+iγ)`), the pole/growth bounds, and the decomposition `F = ∑ m ρ/((σ+iγ)-ρ) + g`
(finite disk zeros `S ∋ β+iγ`, multiplicities `m ρ ≥ 0`, `m(β+iγ) ≥ 1`, `Re g ≥ -B`), the real part of any zero
obeys `β ≤ 1 - 1/(14K)`.  Fully multiplicity-general.  EXIT0 / std-axioms.
-/

open Complex

theorem region_of_zpow_decomp_and_bounds {γ β A B C T₀ T₂ : ℝ} {g F : ℂ} {S : Finset ℂ} {m : ℂ → ℤ}
    (hK : 0 < 3 * A + 4 * B + C) (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hβ : β < 1 + 1 / (2 * (3 * A + 4 * B + C)))
    (hρ₀ : ((β : ℂ) + I * (γ : ℝ)) ∈ S)
    (hleft : ∀ ρ ∈ S, ρ.re < 1 + 1 / (2 * (3 * A + 4 * B + C)))
    (hm : ∀ ρ ∈ S, 0 ≤ m ρ) (hm₀ : 1 ≤ m ((β : ℂ) + I * (γ : ℝ)))
    (hgB : -B ≤ g.re)
    (hdecomp : F = (∑ ρ ∈ S,
      (m ρ : ℂ) / (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ) - ρ)) + g)
    (hdlvp : 0 ≤ 3 * T₀ + 4 * (-F).re + T₂)
    (hpole : T₀ ≤ 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1) + A)
    (hgrowth : T₂ ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hzero : (-F).re ≤ B - 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - β) :=
    re_neg_le_of_zpow_decomp hρ₀ hleft hm hm₀ hβ hgB hdecomp
  exact dlvp_region_of_four_bounds hK hA hB hC hβ hdlvp hpole hzero hgrowth

#print axioms region_of_zpow_decomp_and_bounds
