import Mathlib.Analysis.SpecialFunctions.Gamma.Digamma
import Mathlib.Analysis.SpecialFunctions.Complex.Analytic
import Mathlib.NumberTheory.ArithmeticFunction.VonMangoldt

/-!
# The Weil explicit-formula quadratic form — toward the RH positivity criterion

RH is equivalent to the positivity of the Weil quadratic form
`W(φ) = archimedean(φ) + pole(φ) − prime(φ) = Σ_ρ φ̂(γ_ρ) ≥ 0` (over all suitable test functions).
This module begins a Lean formalization of the *criterion* itself (not just the consequence-side
sign theorem), built on mathlib's `Complex.digamma`.  This is the honest long-term target: a
machine-checkable RH-positivity statement whose components are verified as they are built.

Status (LEDGER §746–750): the RH-hard core is the arithmetic/prime fluctuation; the archimedean
part is the "smooth" piece (Connes' Sonin/prolate domain).  Here we start with the archimedean
weight `Ω(r) = Re[½(ψ(¼+ir/2)+ψ(¼−ir/2))] − log π`, the density against which `|φ̂|²` is integrated.
-/

open Complex

/-- The archimedean weight of the Weil explicit formula (symmetrized real form),
`Ω(r) = Re[½(ψ(¼+ir/2)+ψ(¼−ir/2))] − log π`. -/
noncomputable def archWeight (r : ℝ) : ℝ :=
  (((Complex.digamma (1 / 4 + Complex.I * r / 2)
      + Complex.digamma (1 / 4 - Complex.I * r / 2)) / 2).re) - Real.log Real.pi

/-- The archimedean weight is even in `r` (manifest from the symmetric definition). -/
theorem archWeight_even (r : ℝ) : archWeight (-r) = archWeight r := by
  unfold archWeight
  congr 2
  rw [Complex.ofReal_neg]
  ring_nf


open ArithmeticFunction in
/-- The (finite) prime term of the Weil explicit formula:
`primeₙ(g) = Σ_{n≤N} Λ(n)/√n · g(log n)` — the arithmetic component (LEDGER §749: the RH-hard
prime fluctuation lives here).  `g` is the test-function autocorrelation evaluated at `log n`. -/
noncomputable def weilPrimeTerm (g : ℝ → ℝ) (N : ℕ) : ℝ :=
  ∑ n ∈ Finset.Icc 1 N, (vonMangoldt n / Real.sqrt n) * g (Real.log n)

/-- The prime term vanishes on the zero test function. -/
theorem weilPrimeTerm_zero (N : ℕ) : weilPrimeTerm (fun _ => 0) N = 0 := by
  unfold weilPrimeTerm; simp

/-- The prime term is additive (linear) in the test function. -/
theorem weilPrimeTerm_add (g h : ℝ → ℝ) (N : ℕ) :
    weilPrimeTerm (fun x => g x + h x) N = weilPrimeTerm g N + weilPrimeTerm h N := by
  unfold weilPrimeTerm
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl (fun n _ => ?_)
  ring


/-- The assembled (finite-prime-truncated) Weil functional: given the archimedean integral `A`,
the pole value `P`, and the prime-side test function `g`, `W = A + P − Σ_{n≤N} Λ(n)/√n·g(log n)`.
RH is equivalent to `0 ≤ W` for all test functions (in the `N→∞` limit) — the criterion this
module scaffolds toward.  Here `A`, `P`, `g` are the components a test function induces. -/
noncomputable def weilForm (A P : ℝ) (g : ℝ → ℝ) (N : ℕ) : ℝ := A + P - weilPrimeTerm g N

/-- Weil functional with zero prime data reduces to the archimedean + pole part. -/
theorem weilForm_zero_prime (A P : ℝ) (N : ℕ) : weilForm A P (fun _ => 0) N = A + P := by
  unfold weilForm; rw [weilPrimeTerm_zero]; ring

/-- The Weil functional is nonneg iff the archimedean+pole part dominates the prime term
(the positivity criterion, restated). -/
theorem weilForm_nonneg_iff (A P : ℝ) (g : ℝ → ℝ) (N : ℕ) :
    0 ≤ weilForm A P g N ↔ weilPrimeTerm g N ≤ A + P := by
  unfold weilForm; constructor <;> intro h <;> linarith


/-- The Weil prime term is homogeneous (scalar-linear) in the test function; with
`weilPrimeTerm_add` this gives full linearity — the basis for the Weil form's quadratic structure. -/
theorem weilPrimeTerm_smul (c : ℝ) (g : ℝ → ℝ) (N : ℕ) :
    weilPrimeTerm (fun x => c * g x) N = c * weilPrimeTerm g N := by
  unfold weilPrimeTerm
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl (fun n _ => ?_)
  ring
