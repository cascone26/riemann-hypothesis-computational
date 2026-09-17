import AnalyticAtZeta
import ExplicitFormula
import Mathlib.Analysis.Meromorphic.Order
import Mathlib.Analysis.Meromorphic.Divisor
import Mathlib.NumberTheory.LSeries.Dirichlet

/-!
# §1347 — `divisor ζ U z = 0` for `Re z > 1` (ζ's zeros/poles live in `Re ≤ 1`)

At `z` with `Re z > 1`, ζ is analytic (§1345, `z ≠ 1`) and nonzero (`riemannZeta_ne_zero_of_one_lt_re`), so its
meromorphic order is `0` (`tendsto_ne_zero_iff_meromorphicOrderAt_eq_zero`, ζ tends to a nonzero value), hence its
divisor vanishes (`divisor_apply`).  So the DLVP test point `σ+iγ` (`σ > 1`) is off the divisor support, giving
`d z = 0` (§1346's hypothesis) and — contrapositively — every support point `ρ` has `Re ρ ≤ 1`, so `z ≠ ρ`.
EXIT0 / std-axioms.
-/

open Complex Filter Topology

theorem divisor_riemannZeta_eq_zero_of_one_lt_re {U : Set ℂ} {z : ℂ} (hz : z ∈ U) (hzre : 1 < z.re) :
    MeromorphicOn.divisor riemannZeta U z = 0 := by
  have hmero : MeromorphicOn riemannZeta U := meromorphicOn_riemannZeta.mono_set (Set.subset_univ U)
  have hz1 : z ≠ 1 := by intro h; rw [h] at hzre; simp at hzre
  have hζz : riemannZeta z ≠ 0 := riemannZeta_ne_zero_of_one_lt_re hzre
  have hcont : ContinuousAt riemannZeta z := (analyticAt_riemannZeta hz1).continuousAt
  have hord : meromorphicOrderAt riemannZeta z = 0 :=
    (tendsto_ne_zero_iff_meromorphicOrderAt_eq_zero (hmero z hz)).1
      ⟨riemannZeta z, hζz, hcont.tendsto.mono_left nhdsWithin_le_nhds⟩
  rw [hmero.divisor_apply hz, hord]
  rfl

#print axioms divisor_riemannZeta_eq_zero_of_one_lt_re
