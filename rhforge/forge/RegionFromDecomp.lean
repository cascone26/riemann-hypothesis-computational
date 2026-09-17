import HzeroFromDecomp
import DLVPRegionFromBounds

/-!
# §1326 ★ Region from the ζ decomposition + the three standard bounds (master composition)

Chaining §1314 (partial-fractions decomposition ⟹ `hzero`) into §1311 (four bounds ⟹ region).  At the optimal
abscissa `σ = 1 + 1/(2K)`, `K = 3A+4B+C`, given:
* the `3-4-1` positivity with `T₁ = Re(-F)` (`F = logDeriv ζ (σ+iγ)`),
* the pole bound and the growth bound,
* the partial-fractions decomposition `F = ∑_{ρ∈S} 1/((σ+iγ)-ρ) + g` (finite disk zeros `S ∋ β+iγ`, `Re g ≥ -B`),

the real part of any zero obeys `β ≤ 1 - 1/(14K)`.  Everything is now composed: the whole quantitative zero-free
region reduces to producing that ζ decomposition (via §1316 with `ζ = P·(ζ/P)`) and the standard `T₀,T₂` bounds.
EXIT0 / std-axioms.
-/

open Complex

theorem region_of_decomp_and_bounds {γ β A B C T₀ T₂ : ℝ} {g F : ℂ} {S : Finset ℂ}
    (hK : 0 < 3 * A + 4 * B + C) (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hβ : β < 1 + 1 / (2 * (3 * A + 4 * B + C)))
    (hρ₀ : ((β : ℂ) + I * (γ : ℝ)) ∈ S)
    (hSleft : ∀ ρ ∈ S, ρ.re < 1 + 1 / (2 * (3 * A + 4 * B + C)))
    (hgB : -B ≤ g.re)
    (hdecomp : F = (∑ ρ ∈ S,
      1 / (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ) - ρ)) + g)
    (hdlvp : 0 ≤ 3 * T₀ + 4 * (-F).re + T₂)
    (hpole : T₀ ≤ 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1) + A)
    (hgrowth : T₂ ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hzero : (-F).re ≤ B - 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - β) :=
    re_neg_le_of_decomp hρ₀ hSleft hgB hdecomp
  exact dlvp_region_of_four_bounds hK hA hB hC hβ hdlvp hpole hzero hgrowth

#print axioms region_of_decomp_and_bounds
