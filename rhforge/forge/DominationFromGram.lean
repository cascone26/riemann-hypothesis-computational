import Mathlib
import WeilPositivityDomination

variable {M : Type*} [AddCommGroup M] [Module ℝ M]
variable {V : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]

/-- **Bridge: Gram/pullback structure of the Weil form implies domination (hdom).**
If the Weil bilinear form `A - P` (arch minus prime) is the pullback of an inner
product along a linear map `L` -- `(A - P) v v = inner (L v) (L v)` for all `v` --
then `P v v <= A v v` for all `v`, exactly the forge's `hdom` feeding
`WeilRHFromDomination`.  This is the "RH implies hdom" direction: on RH the explicit
formula writes the Weil form as a genuine Gram form, hence a pullback, hence dominating. -/
theorem domination_from_gram (A P : LinearMap.BilinForm ℝ M) (L : M →ₗ[ℝ] V)
    (h : ∀ v, (A - P) v v = inner ℝ (L v) (L v)) :
    ∀ v, P v v ≤ A v v := by
  intro v
  rw [← weil_positivity_iff_domination A P v, h v]
  exact real_inner_self_nonneg
