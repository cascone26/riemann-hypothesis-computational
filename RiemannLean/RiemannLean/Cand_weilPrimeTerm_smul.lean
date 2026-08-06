import RiemannLean.WeilForm
open ArithmeticFunction
/-- The Weil prime term is homogeneous (scalar-linear) in the test function; with
`weilPrimeTerm_add` this makes it linear — the basis for the Weil form's quadratic structure. -/
theorem weilPrimeTerm_smul (c : ℝ) (g : ℝ → ℝ) (N : ℕ) :
    weilPrimeTerm (fun x => c * g x) N = c * weilPrimeTerm g N := by
  unfold weilPrimeTerm
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl (fun n _ => ?_)
  ring

#print axioms weilPrimeTerm_smul
