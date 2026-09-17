import LFunctionFactorizationAtCenter
import HgcLtOfSphereLFunction
import LFunctionHdecompHgBOfFactor
import BoundLeExplicitLFunction
import GFullNeZero
import Mathlib.NumberTheory.LSeries.Nonvanishing

/-!
# §1508 — `∃g` (half-disk decomposition + `-B ≤ g.re`) for `LFunction χ` from ONLY `g₀`-free data

L-analogue of the ζ CRUX §1424 (cleaner: entire, no `hball`).  Obtain the L factorization on `ball c R'`
internally (§1501), discharge `hgc_lt` from a sphere bound + slack (§1504), get the coordinated decomposition + BC
bound (§1493), and relax the exact `-BOUND ≤ g.re` to `-B ≤ g.re` via the explicit `g₀`-free bound
`BOUND ≤ 4·log(Mζ/‖L c‖)/(R/2)` (§1507) and the scalar `4·log(Mζ/‖L c‖)/(R/2) ≤ B`.  The factor `g₀` never
escapes — the conclusion mentions only `LFunction χ`, its divisor, and the returned `g`.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric Filter Classical

theorem hdecomp_and_hgB_relaxed_of_sphere_LFunction {N : ℕ} [NeZero N]
    {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1) {c : ℂ} {R R' Mζ Ms B : ℝ}
    (hR : 0 < R) (hRR' : R < R') (hcre : 1 < c.re)
    (hζsphere : ∀ z ∈ sphere c R, ‖LFunction χ z‖ ≤ Ms) (hsl : Ms < Mζ) (hMζ : 0 < Mζ)
    (hBexplicit : 4 * Real.log (Mζ / ‖LFunction χ c‖) / (R / 2) ≤ B) :
    ∃ g : ℂ,
      logDeriv (LFunction χ) c
        = (∑ ρ ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
            (fun u => u ∈ closedBall c (R / 2)),
            (MeromorphicOn.divisor (LFunction χ) (ball c R') ρ : ℂ) / (c - ρ)) + g
      ∧ -B ≤ g.re := by
  have hR' : 0 < R' := lt_trans hR hRR'
  obtain ⟨g₀, hga, hgn, hfac⟩ := exists_factorization_LFunction_at_center hχ hR' hcre
  have hζsphereM : ∀ z ∈ sphere c R, ‖LFunction χ z‖ ≤ Mζ := fun z hz => (hζsphere z hz).trans hsl.le
  have hgc_lt := hgc_lt_of_sphere_LFunction hχ hR hRR' hga hfac hζsphere hsl
  obtain ⟨g, hdec, hgBex⟩ :=
    LFunction_hdecomp_and_hgB_of_factor hχ hR hRR' hcre hga hgn hfac hζsphereM hgc_lt
  have hGne : (∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
      (fun u => ¬ u ∈ closedBall c (R / 2)),
      (c - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) * g₀ c ≠ 0 :=
    prod_rest_mul_ne_zero (fun u => (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) g₀
      (fun ρ hρ => (Finset.mem_filter.mp hρ).2) (mem_ball_self (by linarith))
      (hgn ⟨c, mem_ball_self hR'⟩)
  have hGpos : 0 < ‖(∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
      (fun u => ¬ u ∈ closedBall c (R / 2)),
      (c - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) * g₀ c‖ :=
    norm_pos_iff.mpr hGne
  have hZpos : 0 < ‖LFunction χ c‖ :=
    norm_pos_iff.mpr (LFunction_ne_zero_of_one_le_re χ (Or.inl hχ) hcre.le)
  have hBOUND := bound_le_explicit_LFunction hχ hR hRR' hga hfac hMζ hGpos hZpos
  exact ⟨g, hdec, le_trans (neg_le_neg (le_trans hBOUND hBexplicit)) hgBex⟩

#print axioms hdecomp_and_hgB_relaxed_of_sphere_LFunction
