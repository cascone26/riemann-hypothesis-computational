import Mathlib.NumberTheory.LSeries.DirichletContinuation
import Mathlib.NumberTheory.LSeries.Nonvanishing
import Mathlib.Analysis.Meromorphic.Order
import Mathlib.Analysis.Meromorphic.Divisor
import Mathlib.Analysis.Complex.CauchyIntegral

/-!
# §1483 — `divisor (LFunction χ) U z = 0` for `Re z > 1` (`χ ≠ 1`)

L-analogue of ζ §1347.  At `z` with `Re z > 1`, `LFunction χ` is entire (`differentiable_LFunction`) and nonzero
(`DirichletCharacter.LFunction_ne_zero_of_one_le_re`, the Euler-product nonvanishing on `Re ≥ 1`), so its
meromorphic order is `0` (`tendsto_ne_zero_iff_meromorphicOrderAt_eq_zero`), hence its divisor vanishes.  So the
DLVP test point `σ+iγ` (`σ > 1`) is off the divisor support, and — contrapositively — every zero `ρ` has
`Re ρ ≤ 1`.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Filter Topology

theorem divisor_LFunction_eq_zero_of_one_lt_re {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) {U : Set ℂ} {z : ℂ} (hz : z ∈ U) (hzre : 1 < z.re) :
    MeromorphicOn.divisor (LFunction χ) U z = 0 := by
  have hana : AnalyticOnNhd ℂ (LFunction χ) Set.univ :=
    (differentiable_LFunction hχ).differentiableOn.analyticOnNhd isOpen_univ
  have hmero : MeromorphicOn (LFunction χ) U := hana.meromorphicOn.mono_set (Set.subset_univ U)
  have hLz : LFunction χ z ≠ 0 := LFunction_ne_zero_of_one_le_re χ (Or.inl hχ) hzre.le
  have hcont : ContinuousAt (LFunction χ) z := (differentiable_LFunction hχ z).continuousAt
  have hord : meromorphicOrderAt (LFunction χ) z = 0 :=
    (tendsto_ne_zero_iff_meromorphicOrderAt_eq_zero (hmero z hz)).1
      ⟨LFunction χ z, hLz, hcont.tendsto.mono_left nhdsWithin_le_nhds⟩
  rw [hmero.divisor_apply hz, hord]
  rfl

#print axioms divisor_LFunction_eq_zero_of_one_lt_re
