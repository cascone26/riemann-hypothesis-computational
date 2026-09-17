import Mathlib.Analysis.SpecialFunctions.Pow.Complex
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.NumberTheory.DirichletCharacter.Bounds

/-!
# §1511 — twisted truncated-sum modulus bound: `‖∑_{k∈Ioc 1 N} χ(k)·k^{-s}‖ ≤ ∑ k^{-Re s}`

The Dirichlet-character analogue of ζ §1287.  Termwise `‖χ(k)·k^{-s}‖ = ‖χ(k)‖·k^{-Re s} ≤ k^{-Re s}`
(`‖χ(k)‖ ≤ 1`, `DirichletCharacter.norm_le_one`; base `k ≥ 2 > 0`), then the triangle inequality.  This is the
main-sum size estimate for the (forthcoming) truncated Dirichlet-L representation feeding an explicit `L(χ)` growth
bound uniform across `Re s = 1` — the last analytic input for the unconditional L zero-free region.  EXIT0 / std.
-/

open Finset

theorem norm_sum_twist_cpow_neg_le {M : ℕ} [NeZero M] (χ : DirichletCharacter ℂ M) {s : ℂ} {N : ℕ} :
    ‖∑ k ∈ Ioc 1 N, (χ (k : ZMod M)) * (k : ℂ) ^ (-s)‖ ≤ ∑ k ∈ Ioc 1 N, (k : ℝ) ^ (-s.re) := by
  refine le_trans (norm_sum_le _ _) ?_
  refine Finset.sum_le_sum (fun k hk => ?_)
  have hk1 : 1 < k := (Finset.mem_Ioc.mp hk).1
  have hkpos : (0 : ℝ) < (k : ℝ) := by
    exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one (le_of_lt hk1)
  rw [norm_mul, ← Complex.ofReal_natCast k, Complex.norm_cpow_eq_rpow_re_of_pos hkpos, Complex.neg_re]
  calc ‖χ (k : ZMod M)‖ * (k : ℝ) ^ (-s.re)
      ≤ 1 * (k : ℝ) ^ (-s.re) :=
        mul_le_mul_of_nonneg_right (χ.norm_le_one _) (Real.rpow_nonneg hkpos.le _)
    _ = (k : ℝ) ^ (-s.re) := one_mul _

#print axioms norm_sum_twist_cpow_neg_le
