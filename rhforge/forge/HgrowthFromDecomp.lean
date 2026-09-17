import ReNegDecompLe
import ReOneDivSubNonneg

/-!
# §1399 — growth bound `T₂ = Re(-F) ≤ B` from decomposition + BC (`hgrowth` core)

At `z = σ+2iγ` the half-disk decomposition `F = ∑ m/(z-ρ) + g` (`g = logDeriv g_full z`) plus the drop-zeros bound
§1395 (using §1396: `Re(1/(z-ρ)) ≥ 0` since `Re ρ ≤ 1 < σ = Re z`) gives `Re(-F) ≤ -Re g`, and the two-radii BC
bound `-B ≤ Re g` (§1378 at center `σ+2iγ`) gives `-Re g ≤ B`.  Hence `Re(-F) ≤ B` — this is §1336's `hgrowth`
with `C = B ≍ log|t|`, and it is *not* circular (no `1/(σ-1)`).  EXIT0 / std-axioms.
-/

open Complex

theorem hgrowth_from_decomp {S : Finset ℂ} {m : ℂ → ℤ} {z g F : ℂ} {B : ℝ}
    (hF : F = (∑ ρ ∈ S, (m ρ : ℂ) / (z - ρ)) + g) (hm : ∀ ρ ∈ S, 0 ≤ m ρ)
    (hzre : ∀ ρ ∈ S, ρ.re ≤ z.re) (hgB : -B ≤ g.re) :
    (-F).re ≤ B := by
  have h1 : (-F).re ≤ -g.re :=
    re_neg_decomp_le hF hm (fun ρ hρ => re_one_div_sub_nonneg (hzre ρ hρ))
  linarith

#print axioms hgrowth_from_decomp
