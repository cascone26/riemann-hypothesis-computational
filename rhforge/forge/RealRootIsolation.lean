import Mathlib.Algebra.Polynomial.Roots
import Mathlib.Tactic

/-!
# §1808 — Real-root isolation: a real vector selecting one root off a finite set

For `r₀ ∉ rest` (finite multiset of reals) there is a real coefficient vector with `P_c(r₀)=1` and `P_c(s)=0`
for every `s ∈ rest`.  Construction: `P = Q/Q(r₀)` with `Q = ∏_{s∈rest}(X−s)` (real, nonzero at `r₀`).  Feeds
§1807 (`interval_isolated_neg`) to give the full interval converse.  Simpler than the complex isolation
(§1800–1803): no conjugation, just a scalar normalisation.  EXIT0 / std-axioms.
-/

open Polynomial Finset

/-- **Real-root isolation.**  A real vector with `P_c(r₀)=1`, `P_c=0` on `rest`, when `r₀ ∉ rest`. -/
theorem real_root_isolation (rest : Multiset ℝ) (r₀ : ℝ) (h : r₀ ∉ rest) :
    ∃ (n : ℕ) (c : Fin (n + 1) → ℝ),
      (∑ i : Fin (n + 1), c i * r₀ ^ (i : ℕ)) = 1 ∧
      (∀ s ∈ rest, (∑ i : Fin (n + 1), c i * s ^ (i : ℕ)) = 0) := by
  set Q : ℝ[X] := (rest.map (fun s => X - C s)).prod with hQ
  have hQr0 : Q.eval r₀ ≠ 0 := by
    rw [hQ]
    simp only [Polynomial.eval_multiset_prod, Multiset.map_map]
    apply Multiset.prod_ne_zero
    simp only [Multiset.mem_map, Function.comp_apply, Polynomial.eval_sub, Polynomial.eval_X,
      Polynomial.eval_C, not_exists, not_and]
    intro s hs hz; apply h; rw [sub_eq_zero.mp hz]; exact hs
  set P : ℝ[X] := C (Q.eval r₀)⁻¹ * Q with hP
  have hPr0 : P.eval r₀ = 1 := by
    rw [hP, Polynomial.eval_mul, Polynomial.eval_C, inv_mul_cancel₀ hQr0]
  have hPrest : ∀ s ∈ rest, P.eval s = 0 := by
    intro s hs
    rw [hP, Polynomial.eval_mul]
    have hQs : Q.eval s = 0 := by
      rw [hQ]
      simp only [Polynomial.eval_multiset_prod, Multiset.map_map]
      apply Multiset.prod_eq_zero
      rw [Multiset.mem_map]; exact ⟨s, hs, by simp⟩
    rw [hQs, mul_zero]
  refine ⟨P.natDegree, fun i => P.coeff i, ?_, ?_⟩
  · rw [← hPr0, Polynomial.eval_eq_sum_range,
        Fin.sum_univ_eq_sum_range (fun i => P.coeff i * r₀ ^ i) (P.natDegree + 1)]
  · intro s hs
    rw [← hPrest s hs, Polynomial.eval_eq_sum_range,
        Fin.sum_univ_eq_sum_range (fun i => P.coeff i * s ^ i) (P.natDegree + 1)]
