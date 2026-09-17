import HdecompHgBRelaxedOfSphereLFunction
import HgrowthFromDecomp
import DivisorLFunctionNonneg
import DivisorFilterReLeLFunction

/-!
# §1510 — `hgrowth` for `LFunction ψ` from `g₀`-free sphere data: `Re(-logDeriv (LFunction ψ) c) ≤ B`

Combines the CRUX §1508 (half-disk decomposition + BC bound `-B ≤ g.re`, `g₀`-free) with the generic drop-all-zeros
bound §1399 (`hm` from §1479, `hzre` from §1490 with `Re ρ ≤ 1 ≤ Re c`) to give the growth bound
`Re(-logDeriv (LFunction ψ) c) ≤ B` for any nontrivial character `ψ` from sphere data alone.  Applied at `ψ = χ²`
(χ²≠1) and `c = σK+2iγ` this discharges the `hgrowth` hypothesis of the L region §1509.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric

theorem hgrowth_LFunction_from_sphere {N : ℕ} [NeZero N] {ψ : DirichletCharacter ℂ N} (hψ : ψ ≠ 1)
    {c : ℂ} {R R' Mζ Ms B : ℝ} (hR : 0 < R) (hRR' : R < R') (hcre : 1 < c.re)
    (hζsphere : ∀ z ∈ sphere c R, ‖LFunction ψ z‖ ≤ Ms) (hsl : Ms < Mζ) (hMζ : 0 < Mζ)
    (hBexplicit : 4 * Real.log (Mζ / ‖LFunction ψ c‖) / (R / 2) ≤ B) :
    (-logDeriv (LFunction ψ) c).re ≤ B := by
  obtain ⟨g, hdec, hgB⟩ :=
    hdecomp_and_hgB_relaxed_of_sphere_LFunction hψ hR hRR' hcre hζsphere hsl hMζ hBexplicit
  exact hgrowth_from_decomp hdec (fun ρ _ => divisor_LFunction_nonneg_all hψ ρ)
    (fun ρ hρ => divisor_filter_support_re_le_LFunction hψ hcre.le _ hρ) hgB

#print axioms hgrowth_LFunction_from_sphere
