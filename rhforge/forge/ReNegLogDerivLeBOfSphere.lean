import HdecompHgBRelaxedOfSphere
import HgrowthFromDecomp
import DivisorNonnegAll
import DivisorFilterReLt

/-!
# §1429 — `g₀`-free growth bound `Re(-logDeriv ζ c) ≤ B` at any centre `1 < Re c`

Applies the crux §1424 at `c` (obtaining `logDeriv ζ c = ∑ divisor/(c-ρ) + g` with `-B ≤ g.re`, all `g₀`-free),
then drops the half-disk zeros via §1399 (`hgrowth_from_decomp`): each `Re(divisor/(c-ρ)) ≥ 0` since
`Re ρ < 1 < Re c` (§1408), so `Re(-logDeriv ζ c) ≤ -g.re ≤ B`.  Instantiated at `c = σK+2iγ` this supplies
§1426's `hgrowth` with `C := B` (the Borel–Carathéodory bound at the `2iγ` disk, `≍ log|γ|`) — the CORRECT,
fixed-point-compatible growth constant (unlike the pole-order bound §1428).  EXIT0 / std-axioms.
-/

open Complex Metric

theorem re_neg_logDeriv_le_B_of_sphere {c : ℂ} {R R' Mζ Ms B : ℝ} (hR : 0 < R) (hRR' : R < R')
    (hcre : 1 < c.re) (hball : (1 : ℂ) ∉ ball c R')
    (hζsphere : ∀ z ∈ sphere c R, ‖riemannZeta z‖ ≤ Ms) (hsl : Ms < Mζ) (hMζ : 0 < Mζ)
    (hBexplicit : 4 * Real.log (Mζ / ‖riemannZeta c‖) / (R / 2) ≤ B) :
    (-logDeriv riemannZeta c).re ≤ B := by
  obtain ⟨g, hdec, hgB⟩ :=
    hdecomp_and_hgB_relaxed_of_sphere hR hRR' hcre hball hζsphere hsl hMζ hBexplicit
  exact hgrowth_from_decomp hdec (fun ρ _ => divisor_riemannZeta_nonneg_all hball ρ)
    (fun ρ hρ => le_of_lt (divisor_filter_support_re_lt hcre _ hρ)) hgB

#print axioms re_neg_logDeriv_le_B_of_sphere
