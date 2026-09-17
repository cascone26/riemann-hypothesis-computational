import Mathlib.LinearAlgebra.BilinearForm.Basic
import Mathlib.LinearAlgebra.BilinearForm.Properties
import Mathlib.Tactic

/-!
# §1740 — Sector decomposition of a parity-decoupled quadratic form

The formal payoff of the parity block-decoupling (§1733/§1739): once the Weil form's cross term between the
even and odd sectors vanishes, its quadratic form is literally the sum of the two sector quadratic forms.  For a
symmetric bilinear form `B` and a decomposition `v = vₑ + vₒ` with vanishing cross term `B vₑ vₒ = 0`,

`B (vₑ + vₒ) (vₑ + vₒ) = B vₑ vₑ + B vₒ vₒ`.

Hence Weil positivity `B v v ≥ 0` holds iff it holds on each parity sector separately — the RH-relevant
positivity of `arch − prime` decomposes into an even-sector and an odd-sector positivity, which can be analysed
independently.  EXIT0 / std-axioms.
-/

/-- **Quadratic-form sector split.**  A symmetric bilinear form with vanishing cross term is additive across
the two summands. -/
theorem quad_sector_split {M : Type*} [AddCommGroup M] [Module ℝ M]
    (B : LinearMap.BilinForm ℝ M) (hsymm : ∀ x y, B x y = B y x)
    (ve vo : M) (hcross : B ve vo = 0) :
    B (ve + vo) (ve + vo) = B ve ve + B vo vo := by
  simp only [map_add, LinearMap.add_apply]
  rw [hcross, hsymm vo ve, hcross]
  ring
