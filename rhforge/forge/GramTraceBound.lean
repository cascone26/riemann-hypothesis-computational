import CCMWeightGram

/-!
# S2935 — trace (upper) bracket on the zero-side margin (lane 2)

Lower brackets on λ_min came from form bounds (S2929, S2933).  This is the complementary UPPER bracket:
the smallest eigenvalue of a real symmetric matrix is at most its average eigenvalue = trace/n, and the
trace of the weighted Cauchy Gram is explicit:
* `trace_gram` : tr(C·diag(w)·Cᵀ) = Σ_j w_j Σ_k C_kj².
* `exists_eig_mul_card_le_trace` : ∃ i, λ_i · |n| ≤ tr A.
* `ccm_exists_eig_le` : for the CCM weights, ∃ i, λ_i · |n| ≤ (4/L)·‖C‖_F².
Scope: finite truncation, zero side; gives a rigorous ceiling on the margin μ against which the numeric
margin law μ ≈ C·e^{−κλ²} (Report) can be compared.
-/

open Matrix Finset

namespace GramTrace

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq n] [DecidableEq m]

theorem trace_gram (C : Matrix n m ℝ) (w : m → ℝ) :
    (C * diagonal w * Cᵀ).trace = ∑ j, w j * ∑ k, C k j ^ 2 := by
  unfold Matrix.trace
  simp only [Matrix.diag_apply, mul_apply, mul_diagonal, transpose_apply]
  rw [sum_comm]
  refine sum_congr rfl fun j _ => ?_
  rw [mul_sum]
  refine sum_congr rfl fun k _ => ?_
  simp [diagonal_apply]
  ring

theorem exists_eig_mul_card_le_trace [Nonempty n] (A : Matrix n n ℝ) (hA : A.IsHermitian) :
    ∃ i, hA.eigenvalues i * (Fintype.card n : ℝ) ≤ A.trace := by
  by_contra h
  push Not at h
  have htr : A.trace = ∑ i, hA.eigenvalues i := by
    simpa using hA.trace_eq_sum_eigenvalues
  have hlt : ∑ _i : n, A.trace < ∑ i, hA.eigenvalues i * (Fintype.card n : ℝ) :=
    sum_lt_sum_of_nonempty univ_nonempty fun i _ => h i
  rw [sum_const, card_univ, nsmul_eq_mul, ← sum_mul, ← htr] at hlt
  linarith [mul_comm (Fintype.card n : ℝ) A.trace]

theorem ccm_exists_eig_le [Nonempty n] {L : ℝ} (hL : 0 < L) (C : Matrix n m ℝ) (γ : m → ℝ) :
    ∃ i, (CCMWeight.ccm_gram_hermitian hL C γ).eigenvalues i * (Fintype.card n : ℝ)
      ≤ 4 / L * ∑ j, ∑ k, C k j ^ 2 := by
  obtain ⟨i, hi⟩ := exists_eig_mul_card_le_trace _ (CCMWeight.ccm_gram_hermitian hL C γ)
  refine ⟨i, hi.trans ?_⟩
  rw [trace_gram, mul_sum]
  exact sum_le_sum fun j _ => mul_le_mul_of_nonneg_right (CCMWeight.ccmW_le hL (γ j))
    (sum_nonneg fun k _ => sq_nonneg _)

end GramTrace

#print axioms GramTrace.trace_gram
#print axioms GramTrace.exists_eig_mul_card_le_trace
#print axioms GramTrace.ccm_exists_eig_le
