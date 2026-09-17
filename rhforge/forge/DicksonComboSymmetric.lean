import DicksonSum

/-!
# §1784 — Dickson combination is a symmetric (palindromic) Laurent polynomial

The reduction heart of the Chebyshev on-circle route.  §1783 gave `D_n(z+z⁻¹) = zⁿ+z⁻ⁿ`.  Here we assemble
a general combination: for coefficients `c₀,…,c_m` and `z ≠ 0`,

    zᵐ · Σ_{j≤m} cⱼ · D_j(z+z⁻¹)  =  Σ_{j≤m} cⱼ · (z^{m+j} + z^{m-j}).

The right-hand side ranges over powers `z^{m-j} … z^{m+j}` symmetrically about `zᵐ`, i.e. it is exactly a
palindromic polynomial `P(z) = zᵐ·Q(z+z⁻¹)` with `Q = Σ cⱼ D_j`.  This is the identity `P(z) = zᵐ Q(z+z⁻¹)`
underlying the degree-reduction: with §1781/§1782 it turns "P's zeros on `|z|=1`" into "Q's roots real in
`[-2,2]`".  EXIT0 / std-axioms.
-/

open Polynomial

variable {K : Type*} [Field K]

/-- **Dickson combination = symmetric Laurent polynomial.**  For `z ≠ 0`, multiplying a Dickson combination
`Σ cⱼ D_j(z+z⁻¹)` by `zᵐ` yields the palindromic power sum `Σ cⱼ (z^{m+j} + z^{m-j})`. -/
theorem dickson_combo_symmetric (c : ℕ → K) (z : K) (hz : z ≠ 0) (m : ℕ) :
    z ^ m * (∑ j ∈ Finset.range (m + 1), c j * (dickson j).eval (z + z⁻¹))
      = ∑ j ∈ Finset.range (m + 1), c j * (z ^ (m + j) + z ^ (m - j)) := by
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro j hj
  have hjm : j ≤ m := Nat.lt_succ_iff.mp (Finset.mem_range.mp hj)
  have hzj : z ^ j ≠ 0 := pow_ne_zero j hz
  rw [dickson_eval_sum z hz j, inv_pow, pow_add]
  have e2 : z ^ (m - j) = z ^ m * (z ^ j)⁻¹ := by
    field_simp
    rw [← pow_add, Nat.sub_add_cancel hjm]
  rw [e2]
  ring
