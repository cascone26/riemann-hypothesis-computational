import P2KLPoly

open scoped MatrixOrder Matrix

set_option maxHeartbeats 1000000

/-!
# Off-critical-line lane: `P'_{2K,L}(ν_{m,L}) = ∏_{j≠m}(ν_{m,L}−ν_{j,L})` (feeds Eq 7.13)

Shi Eq (7.13)'s interpolation coefficients `c^{off}_{n,2K,L} := Q^{off}_K(ν_{n,L})/P'_{2K,L}(ν_{n,L})`
need the derivative of §2366-2367's `P_2KL` at each frequency node. This brick computes it exactly:
via `HasDerivAt.finsetProd` (differentiating the finite product termwise) followed by `Finset.
sum_eq_single` (every term of the resulting sum except the one indexed by `m` itself contains the
vanishing factor `(ν_{m,L}−ν_{m,L})=0`, since every OTHER term's product ranges over `univ.erase n`
for `n≠m`, which still contains `m`). Search-before-build: grepped for `P_2KL_hasDerivAt`/
`P_2KL_deriv` — fresh. Fully general, no zeta/RH dependence. -/
theorem P_2KL_hasDerivAt_node (K : ℕ) (L : ℝ) (m : Fin (2 * (2 * K) + 1)) :
    HasDerivAt (P_2KL K L)
      (∏ j ∈ Finset.univ.erase m,
        ((2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L)
          - (2 * Real.pi * (((j : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L)))
      (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L) := by
  unfold P_2KL
  set x : ℝ := 2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L with hxdef
  have hfactor : ∀ n : Fin (2 * (2 * K) + 1),
      HasDerivAt (fun w : ℝ => w - 2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L) 1 x :=
    fun n => (hasDerivAt_id x).sub_const _
  have hprod := HasDerivAt.finsetProd (u := (Finset.univ : Finset (Fin (2 * (2 * K) + 1))))
    (fun n _ => hfactor n)
  have hsimp : (∑ n : Fin (2 * (2 * K) + 1),
      (∏ j ∈ Finset.univ.erase n, (x - 2 * Real.pi * (((j : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L)) • (1:ℝ))
      = ∏ j ∈ Finset.univ.erase m, (x - 2 * Real.pi * (((j : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L) := by
    rw [Finset.sum_eq_single m]
    · simp
    · intro n _ hnm
      have hmem : m ∈ Finset.univ.erase n := Finset.mem_erase.mpr ⟨Ne.symm hnm, Finset.mem_univ m⟩
      have hzero : x - 2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L = 0 := by
        rw [hxdef]; ring
      rw [Finset.prod_eq_zero hmem hzero]
      simp
    · intro h; exact absurd (Finset.mem_univ m) h
  rw [hsimp] at hprod
  have heq : (∏ i : Fin (2 * (2 * K) + 1),
      fun w : ℝ => w - 2 * Real.pi * (((i : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L)
      = fun w => ∏ n : Fin (2 * (2 * K) + 1), (w - 2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L) := by
    funext w
    rw [Finset.prod_apply]
  rw [heq] at hprod
  exact hprod

#print axioms P_2KL_hasDerivAt_node
