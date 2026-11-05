import FreqNodeReflect

open scoped MatrixOrder Matrix

/-!
# Off-critical-line lane: `P_{2K,L}` is odd under negation (step 3 of the matrix-connection plan)

Search-before-build: grepped `P_2KL_neg`/`P_2KL_odd`/`P2KL.*neg` — fresh. The frequency-node set
`{ν_{n,L}}` is symmetric under negation (§2392-2394's `freqReflect`/`freqNode_reflect_eq_neg`), and
the index set has odd cardinality `4K+1`, so `P_{2K,L}(-w) = (-1)^{4K+1}·∏(w+ν_n) =
-∏(w-ν_{reflect(n)}) = -P_{2K,L}(w)` after reindexing by the reflection bijection. This is the key
fact needed (together with `Q_off_K_even`) to show `c_off` is invariant under `freqReflect`,
completing step 3 of §2391's 4-step plan. Fully general, no zeta/RH dependence. -/
theorem P_2KL_odd (K : ℕ) (L : ℝ) (w : ℝ) :
    P_2KL K L (-w) = -P_2KL K L w := by
  unfold P_2KL
  have step1 : ∀ n : Fin (2 * (2 * K) + 1),
      -w - 2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L
        = (-1) * (w + 2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L) := by
    intro n; ring
  simp_rw [step1]
  rw [Finset.prod_mul_distrib, Finset.prod_const, Finset.card_univ, Fintype.card_fin]
  have hodd : Odd (2 * (2 * K) + 1) := ⟨2 * K, rfl⟩
  rw [hodd.neg_one_pow]
  have step2 : ∏ n : Fin (2 * (2 * K) + 1),
      (w + 2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L)
        = ∏ n : Fin (2 * (2 * K) + 1),
      (w - 2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L) := by
    rw [← Function.Bijective.prod_comp (freqReflect_bijective K)
      (fun n => w - 2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L)]
    apply Finset.prod_congr rfl
    intro n _
    have h := freqNode_reflect_eq_neg K L n
    linarith [h]
  rw [step2]
  ring

#print axioms P_2KL_odd
