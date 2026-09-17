import Mathlib.Analysis.SpecialFunctions.Pow.Complex
import Mathlib.Analysis.SpecialFunctions.Pow.Real

/-!
# §1287 — Modulus bound on the truncated main sum: `‖∑_{k≤N} k^{-s}‖ ≤ ∑_{k≤N} k^{-Re s}`

First size estimate on the finite Euler–Maclaurin main sum (§1286).  Termwise `‖k^{-s}‖ = k^{-Re s}` for
`k ≥ 1` (base positive real), then the triangle inequality.  With `Re s = σ` near `1` the RHS `∑_{k≤N} k^{-σ}`
is `≍ log N`; taking `N ≍ |t|` this is the `log|t|` growth input to §1269's conditional bound.
EXIT0 / std-axioms.
-/

open Finset

theorem norm_sum_cpow_neg_le {N : ℕ} {s : ℂ} :
    ‖∑ k ∈ Ioc 1 N, (k : ℂ) ^ (-s)‖ ≤ ∑ k ∈ Ioc 1 N, (k : ℝ) ^ (-s.re) := by
  refine le_trans (norm_sum_le _ _) ?_
  refine Finset.sum_le_sum (fun k hk => ?_)
  have hk1 : 1 < k := (Finset.mem_Ioc.mp hk).1
  have hkpos : (0 : ℝ) < (k : ℝ) := by exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one (le_of_lt hk1)
  rw [← Complex.ofReal_natCast k, Complex.norm_cpow_eq_rpow_re_of_pos hkpos, Complex.neg_re]

#print axioms norm_sum_cpow_neg_le
