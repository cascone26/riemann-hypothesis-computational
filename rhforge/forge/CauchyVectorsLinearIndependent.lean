import Mathlib

open Polynomial

theorem cauchy_vectors_linear_independent (m : ℕ) (hm : 1 ≤ m) (a : Fin m → ℂ)
    (ν : Fin (m + 1) → ℝ)
    (ha_inj : Function.Injective a) (ha_im : ∀ j, (a j).im ≠ 0)
    (hν_inj : Function.Injective ν)
    (c : Fin m → ℂ) (hc : ∀ n : Fin (m + 1), ∑ j, c j / (a j + (ν n : ℂ)) = 0) :
    ∀ j, c j = 0 := by
  classical
  set P : ℂ[X] := ∑ j : Fin m, C (c j) * ∏ i ∈ Finset.univ.erase j, (X - C (-(a i))) with hPdef
  have hne : ∀ (n : Fin (m + 1)) (j : Fin m), a j + (ν n : ℂ) ≠ 0 := by
    intro n j heq
    apply ha_im j
    have : (a j + (ν n : ℂ)).im = 0 := by rw [heq]; simp
    simpa using this
  have hPeval : ∀ n : Fin (m + 1), P.eval (ν n : ℂ) = 0 := by
    intro n
    have hkey : ∀ j : Fin m,
        (C (c j) * ∏ i ∈ Finset.univ.erase j, (X - C (-(a i)))).eval (ν n : ℂ)
        = c j * ∏ i ∈ Finset.univ.erase j, (a i + (ν n : ℂ)) := by
      intro j
      simp only [eval_mul, eval_C, eval_prod, eval_sub, eval_X]
      congr 1
      apply Finset.prod_congr rfl
      intro i _
      ring
    have hsum_eval : P.eval (ν n : ℂ)
        = ∑ j : Fin m, c j * ∏ i ∈ Finset.univ.erase j, (a i + (ν n : ℂ)) := by
      rw [hPdef, eval_finsetSum]
      exact Finset.sum_congr rfl (fun j _ => hkey j)
    rw [hsum_eval]
    have hmul : (∑ j, c j / (a j + (ν n : ℂ))) * ∏ i, (a i + (ν n : ℂ)) = 0 := by
      rw [hc n]; ring
    have hexpand : (∑ j, c j / (a j + (ν n : ℂ))) * ∏ i, (a i + (ν n : ℂ))
        = ∑ j, c j * ∏ i ∈ Finset.univ.erase j, (a i + (ν n : ℂ)) := by
      rw [Finset.sum_mul]
      apply Finset.sum_congr rfl
      intro j _
      rw [← Finset.mul_prod_erase Finset.univ (fun i => a i + (ν n : ℂ)) (Finset.mem_univ j),
        ← mul_assoc, div_mul_cancel₀ _ (hne n j)]
    rw [← hexpand]
    exact hmul
  have hdeg : P.natDegree ≤ m - 1 := by
    rw [hPdef]
    apply Polynomial.natDegree_sum_le_of_forall_le
    intro j _
    calc (C (c j) * ∏ i ∈ Finset.univ.erase j, (X - C (-(a i)))).natDegree
        ≤ (C (c j)).natDegree + (∏ i ∈ Finset.univ.erase j, (X - C (-(a i)))).natDegree :=
          Polynomial.natDegree_mul_le
      _ = (∏ i ∈ Finset.univ.erase j, (X - C (-(a i)))).natDegree := by
          simp [Polynomial.natDegree_C]
      _ = (Finset.univ.erase j).card := by
          rw [Polynomial.natDegree_prod_of_monic _ _ (fun i _ => Polynomial.monic_X_sub_C (-(a i)))]
          simp [Polynomial.natDegree_X_sub_C]
      _ = m - 1 := by
          rw [Finset.card_erase_of_mem (Finset.mem_univ j), Finset.card_univ, Fintype.card_fin]
  have hs_card : (Finset.image (fun n : Fin (m + 1) => (ν n : ℂ)) Finset.univ).card = m + 1 := by
    rw [Finset.card_image_of_injective _ (fun x y h => hν_inj (by exact_mod_cast h))]
    simp
  have hPzero : P = 0 := by
    apply Polynomial.eq_zero_of_degree_lt_of_eval_finset_eq_zero
      (Finset.image (fun n : Fin (m + 1) => (ν n : ℂ)) Finset.univ)
    · rw [hs_card]
      calc P.degree ≤ (P.natDegree : WithBot ℕ) := Polynomial.degree_le_natDegree
        _ ≤ ((m - 1 : ℕ) : WithBot ℕ) := by exact_mod_cast hdeg
        _ < ((m + 1 : ℕ) : WithBot ℕ) := by
            have : m - 1 < m + 1 := by omega
            exact_mod_cast this
    · intro x hx
      simp only [Finset.mem_image, Finset.mem_univ, true_and] at hx
      obtain ⟨n, hn⟩ := hx
      rw [← hn]
      exact hPeval n
  intro j
  have hPj : P.eval (-(a j)) = 0 := by rw [hPzero]; simp
  have hterms : P.eval (-(a j))
      = ∑ i : Fin m, c i * ∏ k ∈ Finset.univ.erase i, (a k + (-(a j))) := by
    rw [hPdef, eval_finsetSum]
    apply Finset.sum_congr rfl
    intro i _
    simp only [eval_mul, eval_C, eval_prod, eval_sub, eval_X]
    congr 1
    apply Finset.prod_congr rfl
    intro k _
    ring
  rw [hterms] at hPj
  have hvanish : ∀ i : Fin m, i ≠ j →
      c i * ∏ k ∈ Finset.univ.erase i, (a k + (-(a j))) = 0 := by
    intro i hij
    have hmem : j ∈ Finset.univ.erase i := Finset.mem_erase.mpr ⟨Ne.symm hij, Finset.mem_univ j⟩
    have hz : (∏ k ∈ Finset.univ.erase i, (a k + (-(a j)))) = 0 := by
      apply Finset.prod_eq_zero hmem
      ring
    rw [hz]; ring
  rw [Finset.sum_eq_single j (fun i _ hij => hvanish i hij) (by simp)] at hPj
  have hprodne : (∏ k ∈ Finset.univ.erase j, (a k + (-(a j)))) ≠ 0 := by
    apply Finset.prod_ne_zero_iff.mpr
    intro k hk
    simp only [Finset.mem_erase] at hk
    have : a k ≠ a j := fun h => hk.1 (ha_inj h)
    intro hzero
    apply this
    linear_combination hzero
  exact (mul_eq_zero.mp hPj).resolve_right hprodne

#print axioms cauchy_vectors_linear_independent
