import P2KLOdd

open scoped MatrixOrder Matrix

set_option maxHeartbeats 1000000

/-!
# Off-critical-line lane: `P'_{2K,L}` is even (derivative of an odd function is even)

Search-before-build: grepped `deriv.*odd`/`odd.*deriv`/`P_2KL_deriv_even`/`P_2KL_hasDerivAt_neg` —
fresh. Step 3-c of §2391's 4-step plan. First a fully GENERAL, non-RH-specific lemma: if
`F(-x)=-F(x)` for all `x` and `F` is differentiable everywhere with derivative function `F'`, then
`F'` is even (`F'(-x)=F'(x)`) — via the chain rule on `F∘Neg.neg` (two ways to compute its
derivative: directly, and via `F(-x)=-F(x)` rewritten as `-F`) plus `HasDerivAt.unique`. Then a
GENERAL (any real point, not just nodes) `HasDerivAt` fact for `P_2KL` via the same
`HasDerivAt.finsetProd` Leibniz-rule technique as §2368's node-specific version (simpler here since
there is no erase-collapsing needed). Combining the two with §2395's `P_2KL_odd` gives
`P_2KL`'s derivative is even everywhere — the key remaining ingredient for
`c_off_freqReflect_invariant`. Fully general, no zeta/RH dependence. -/
theorem hasDerivAt_odd_deriv_even (F : ℝ → ℝ) (F' : ℝ → ℝ)
    (hF : ∀ x, HasDerivAt F (F' x) x) (hodd : ∀ x, F (-x) = -F x) (x : ℝ) :
    F' (-x) = F' x := by
  have hcomp0 := (hF (-x)).comp x (hasDerivAt_neg x)
  have heq : (F ∘ Neg.neg) = fun w => -F w := by
    funext w
    simp only [Function.comp_apply]
    exact hodd w
  rw [heq] at hcomp0
  have hneg : HasDerivAt (fun w => -F w) (-F' x) x := (hF x).neg
  have hu := hcomp0.unique hneg
  linarith [hu]

theorem P_2KL_hasDerivAt (K : ℕ) (L : ℝ) (w : ℝ) :
    HasDerivAt (P_2KL K L)
      (∑ n : Fin (2 * (2 * K) + 1), ∏ j ∈ Finset.univ.erase n,
        (w - 2 * Real.pi * (((j : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L))
      w := by
  unfold P_2KL
  have hfactor : ∀ n : Fin (2 * (2 * K) + 1),
      HasDerivAt (fun v : ℝ => v - 2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L) 1 w :=
    fun n => (hasDerivAt_id w).sub_const _
  have hprod := HasDerivAt.finsetProd (u := (Finset.univ : Finset (Fin (2 * (2 * K) + 1))))
    (fun n _ => hfactor n)
  simp only [smul_eq_mul, mul_one] at hprod
  have heq : (∏ i : Fin (2 * (2 * K) + 1),
      fun v : ℝ => v - 2 * Real.pi * (((i : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L)
      = fun v => ∏ n : Fin (2 * (2 * K) + 1),
        (v - 2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L) := by
    funext v
    rw [Finset.prod_apply]
  rw [heq] at hprod
  exact hprod

theorem P_2KL_deriv_even (K : ℕ) (L : ℝ) (w : ℝ) :
    deriv (P_2KL K L) (-w) = deriv (P_2KL K L) w := by
  have hF : ∀ x, HasDerivAt (P_2KL K L) (deriv (P_2KL K L) x) x :=
    fun x => (P_2KL_hasDerivAt K L x).deriv ▸ (P_2KL_hasDerivAt K L x)
  exact hasDerivAt_odd_deriv_even (P_2KL K L) (deriv (P_2KL K L)) hF (P_2KL_odd K L) w

#print axioms hasDerivAt_odd_deriv_even
#print axioms P_2KL_hasDerivAt
#print axioms P_2KL_deriv_even
