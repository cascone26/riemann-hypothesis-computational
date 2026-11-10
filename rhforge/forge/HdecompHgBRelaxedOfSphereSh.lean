import ZetaFactorizationAtCenter
import HgcLtOfSphere
import ZetaHdecompHgBOfFactorSh
import BoundLeExplicitSh
import GFullNeZero
import Mathlib.NumberTheory.LSeries.Dirichlet

/-! [S3020 SHARP COPY: Borel–Carathéodory constant 4 → 2 via SharpBC (S3018); names suffixed _sh]

# §1424 — `∃g` (half-disk decomposition + `-B ≤ g.re`) from ONLY `g₀`-free data

The crux: obtain the ζ factorization on `ball c R'` internally (§1423), discharge `hgc_lt` from a sphere bound +
slack (§1414), get the coordinated decomposition + exact BC bound (§1410), and relax the exact `-BOUND ≤ g.re` to
`-B ≤ g.re` using the explicit `g₀`-free bound `BOUND ≤ 4·log(Mζ/‖ζ c‖)/(R/2)` (§1419) and the scalar hypothesis
`4·log(Mζ/‖ζ c‖)/(R/2) ≤ B`.  The factor `g₀` never escapes: the conclusion mentions only ζ, the divisor, and the
returned factor value `g`.  This is exactly §1409's `hdecomp`+`hgB` pair, now with **no reference to the
factorization**.  EXIT0 / std-axioms.
-/

open Complex Metric Filter Classical

theorem hdecomp_and_hgB_relaxed_of_sphere_sh {c : ℂ} {R R' Mζ Ms B : ℝ} (hR : 0 < R) (hRR' : R < R')
    (hcre : 1 < c.re) (hball : (1 : ℂ) ∉ ball c R')
    (hζsphere : ∀ z ∈ sphere c R, ‖riemannZeta z‖ ≤ Ms) (hsl : Ms < Mζ) (hMζ : 0 < Mζ)
    (hBexplicit : 2 * Real.log (Mζ / ‖riemannZeta c‖) / (R / 2) ≤ B) :
    ∃ g : ℂ,
      logDeriv riemannZeta c
        = (∑ ρ ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
            (fun u => u ∈ closedBall c (R / 2)),
            (MeromorphicOn.divisor riemannZeta (ball c R') ρ : ℂ) / (c - ρ)) + g
      ∧ -B ≤ g.re := by
  have hR' : 0 < R' := lt_trans hR hRR'
  obtain ⟨g₀, hga, hgn, hfac⟩ := exists_factorization_riemannZeta_at_center hR' hcre
  have hζsphereM : ∀ z ∈ sphere c R, ‖riemannZeta z‖ ≤ Mζ := fun z hz => (hζsphere z hz).trans hsl.le
  have hgc_lt := hgc_lt_of_sphere hR hRR' hball hga hfac hζsphere hsl
  obtain ⟨g, hdec, hgBex⟩ :=
    zeta_hdecomp_and_hgB_of_factor_sh hR hRR' hcre hball hga hgn hfac hζsphereM hgc_lt
  have hGne : (∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
      (fun u => ¬ u ∈ closedBall c (R / 2)),
      (c - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) * g₀ c ≠ 0 :=
    prod_rest_mul_ne_zero (fun u => (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) g₀
      (fun ρ hρ => (Finset.mem_filter.mp hρ).2) (mem_ball_self (by linarith))
      (hgn ⟨c, mem_ball_self hR'⟩)
  have hGpos : 0 < ‖(∏ u ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter
      (fun u => ¬ u ∈ closedBall c (R / 2)),
      (c - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R') u).toNat) * g₀ c‖ :=
    norm_pos_iff.mpr hGne
  have hZpos : 0 < ‖riemannZeta c‖ := norm_pos_iff.mpr (riemannZeta_ne_zero_of_one_lt_re hcre)
  have hBOUND := bound_le_explicit_sh hR hRR' hball hga hfac hMζ hGpos hZpos
  exact ⟨g, hdec, le_trans (neg_le_neg (le_trans hBOUND hBexplicit)) hgBex⟩

#print axioms hdecomp_and_hgB_relaxed_of_sphere_sh
