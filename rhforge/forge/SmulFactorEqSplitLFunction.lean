import SmulSplitEval
import DivisorLFunctionNonneg
import LFunctionDivisorFinite

/-!
# §1480 — `LFunction χ`'s factorized-rational factor equals `P_half · g_full` (as functions)

L-analogue of ζ §1388.  Specialising the generic split-evaluation `smul_split_eval_nat` (§1386) to `LFunction χ`
— with finite divisor support §1475 and divisor `≥ 0` §1479 — and taking `funext`: the factorized-rational smul
from §1476 equals, *as a function*, `fun s => P_half(s) · g_full(s)`, where `P_half`/`g_full` split the divisor
support by the predicate `p`.  This rewrites §1476's codiscrete factorization into `P_half · g_full` form before
the pointwise (sphere) upgrade.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric

open Classical in
theorem smul_factor_eq_split_LFunction {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {c : ℂ} {R : ℝ} {g : ℂ → ℂ} (p : ℂ → Prop) :
    ((∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R) u)) • g)
      = fun s =>
        (∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R).toFinset.filter p,
            (s - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R) u).toNat)
          * ((∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R).toFinset.filter (fun u => ¬ p u),
              (s - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R) u).toNat) * g s) := by
  funext s
  exact smul_split_eval_nat (divisor_LFunction_ball_support_finite hχ c R)
    (divisor_LFunction_nonneg_all hχ) s p

#print axioms smul_factor_eq_split_LFunction
