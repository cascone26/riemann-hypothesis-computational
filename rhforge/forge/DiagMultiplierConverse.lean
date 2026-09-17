import Mathlib.Tactic

/-!
# §1762 — Converse engine of the RH criteria: a negative diagonal entry breaks positivity

The RH-direction bricks (§1759 `diagForm_nonneg`, §1761 capstone) show `M ≥ 0` (real zeros, positive
multiplicities) ⟹ the Weil form is PSD.  The CONVERSE direction of Suzuki's criteria (¬RH ⟹ ∃ψ with
`⟨ψ,ψ⟩_W < 0`) rests on the dual algebraic fact: in the zero-basis the Weil form is diagonal,
`W(ψ,ψ) = Σ_γ m_γ |c_γ(ψ)|²`, and an off-line zero `γ₀` contributes a NEGATIVE weight `m_{γ₀} < 0`
(Suzuki Thm 4.4 converse: the `−m_{γ₀}` term).  A diagonal form with any negative entry fails to be PSD —
concentrate the input on that coordinate and the form goes negative.

`M i < 0  ⟹  diagForm M (𝐞_i) = M i < 0`  — so the form is not positive semidefinite.

This is the exact abstract engine of the converse: a single off-line zero (negative diagonal weight) produces
a negative Weil-form value, contradicting positivity.  Companion to §1759 (`M ≥ 0 ⟹ PSD`); together they give
the diagonal-form both-directions dictionary underlying RH ⟺ (all weights `≥ 0`).  EXIT0 / std-axioms.
-/

open Finset

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The diagonal multiplier form `q(x) = ∑ i, M i · (x i)²` (as in §1759). -/
def diagForm (M x : ι → ℝ) : ℝ := ∑ i, M i * (x i * x i)

/-- **A negative diagonal entry is exposed by its own basis vector.** -/
theorem diagForm_single (M : ι → ℝ) (i : ι) :
    diagForm M (Pi.single i 1) = M i := by
  unfold diagForm
  rw [Finset.sum_eq_single i]
  · rw [Pi.single_eq_same]; ring
  · intro j _ hj; rw [Pi.single_eq_of_ne hj]; ring
  · intro h; exact absurd (Finset.mem_univ i) h

/-- **Converse engine.**  If any diagonal weight is negative, the form is not PSD: an explicit input
(the basis vector at that coordinate) makes it strictly negative. -/
theorem diagForm_neg_of_neg_entry (M : ι → ℝ) (i : ι) (hi : M i < 0) :
    diagForm M (Pi.single i 1) < 0 := by
  rw [diagForm_single]; exact hi

/-- **Not-PSD from a negative entry** (packaged existential — the shape used by the RH converse). -/
theorem diagForm_not_psd_of_neg_entry (M : ι → ℝ) (i : ι) (hi : M i < 0) :
    ∃ x, diagForm M x < 0 :=
  ⟨Pi.single i 1, diagForm_neg_of_neg_entry M i hi⟩
