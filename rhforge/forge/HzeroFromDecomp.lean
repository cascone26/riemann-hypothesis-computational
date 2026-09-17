import ReNegPartialFractions
import ReOneDivSameIm

/-!
# §1314 ★ `hzero` from the Hadamard decomposition — the exact §1311 input, reduced to ONE analytic fact

Combining §1312 (drop-zeros → `Re(-(∑+g)) ≤ B - (1/(s-ρ₀)).re`) and §1313 (`(1/(s-ρ₀)).re = 1/(σ-β)` at same
height), the DLVP `hzero` input is obtained from the sole open analytic fact — the partial-fractions / Hadamard
decomposition `F = (∑_{ρ∈S} 1/(s-ρ)) + g` of `-ζ'/ζ` on a disk (finite zero set `S ∋ ρ₀=β+iγ`, bounded
remainder `Re g ≥ -B`).  For any `F` admitting that decomposition at `s = σ+iγ`:

  `Re(-F) ≤ B - 1/(σ-β)`.

Feeding `F = ζ'/ζ(σ+iγ)` (i.e. `T₁ = Re(-ζ'/ζ) ≤ B - 1/(σ-β)`) into §1311 closes the quantitative region.
The ONLY thing still owed is the decomposition itself (Borel–Carathéodory local expansion).  EXIT0 / std-axioms.
-/

open Complex

theorem re_neg_le_of_decomp {σ β γ B : ℝ} {g F : ℂ} {S : Finset ℂ}
    (hρ₀ : ((β : ℂ) + I * (γ : ℝ)) ∈ S)
    (hSleft : ∀ ρ ∈ S, ρ.re < σ) (hgB : -B ≤ g.re)
    (hdecomp : F = (∑ ρ ∈ S, 1 / (((σ : ℂ) + I * (γ : ℝ)) - ρ)) + g) :
    (-F).re ≤ B - 1 / (σ - β) := by
  have hsre : (((σ : ℂ) + I * (γ : ℝ))).re = σ := by simp
  have hleft' : ∀ ρ ∈ S, ρ.re < (((σ : ℂ) + I * (γ : ℝ))).re := by
    intro ρ hρ; rw [hsre]; exact hSleft ρ hρ
  have h12 := re_neg_partial_fractions_le (s := (σ : ℂ) + I * (γ : ℝ))
    (ρ₀ := (β : ℂ) + I * (γ : ℝ)) hρ₀ hleft' hgB
  rw [re_one_div_sub_same_im σ β γ] at h12
  rw [hdecomp]
  exact h12

#print axioms re_neg_le_of_decomp
