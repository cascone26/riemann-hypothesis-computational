import QOffKPoly

open scoped MatrixOrder Matrix

/-!
# Off-critical-line lane: `P_{2K,L}(w)` (Eq 7.10), the frequency-node interpolation polynomial

Shi Eq (7.10): `P_{2K,L}(w) := ∏_{n=-2K}^{2K} (w−ν_{n,L})` — the companion polynomial to §2364-2365's
`Q^{off}_K`, whose ratio `Q^{off}_K(ν_{n,L})/P'_{2K,L}(ν_{n,L})` (Eq 7.13) gives the interpolation
coefficients `c^{off}_{n,2K,L}`. This brick defines `P_2KL` (using this session's established
`ν_{n,L}=2π(n-N)/L` convention with `N:=2K`, matching every other object built this session over the
same frequency-node index) and proves the essential root fact: `P_2KL` vanishes at each frequency
node `ν_{n,L}`. Search-before-build: grepped for `P_2K_L`/`P2KL`/`freqNodeProduct` — fresh (only this
file's own forward-reference comments existed). Fully general, no zeta/RH dependence. -/
noncomputable def P_2KL (K : ℕ) (L : ℝ) (w : ℝ) : ℝ :=
  ∏ n : Fin (2 * (2 * K) + 1), (w - 2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L)

theorem P_2KL_eq_zero_at_node (K : ℕ) (L : ℝ) (m : Fin (2 * (2 * K) + 1)) :
    P_2KL K L (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L) = 0 := by
  unfold P_2KL
  apply Finset.prod_eq_zero (Finset.mem_univ m)
  ring

#print axioms P_2KL_eq_zero_at_node
