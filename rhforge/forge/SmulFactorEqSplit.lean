import SmulSplitEval
import DivisorNonnegAll
import ZetaDivisorFinite

/-!
# §1388 — ζ's factorized-rational factor equals `P_half · g_full` (as functions)

Specialising §1386 to ζ (with `divisor ≥ 0` from §1387) and taking `funext`: the factorized-rational smul from
§1343 equals, *as a function*, `fun s => P_half(s) · g_full(s)`, where `P_half`/`g_full` split the divisor support
by the predicate `p`.  This is the function-level identity that lets §1343's codiscrete factorization be rewritten
into `P_half · g_full` form before the pointwise (sphere) upgrade §1385.  EXIT0 / std-axioms.
-/

open Metric

open Classical in
theorem smul_factor_eq_split {c : ℂ} {R : ℝ} (hball : (1 : ℂ) ∉ ball c R) {g : ℂ → ℂ}
    (p : ℂ → Prop) :
    ((∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R) u)) • g)
      = fun s =>
        (∏ u ∈ (divisor_riemannZeta_ball_support_finite c R).toFinset.filter p,
            (s - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R) u).toNat)
          * ((∏ u ∈ (divisor_riemannZeta_ball_support_finite c R).toFinset.filter (fun u => ¬ p u),
              (s - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R) u).toNat) * g s) := by
  funext s
  exact smul_split_eval_nat (divisor_riemannZeta_ball_support_finite c R)
    (divisor_riemannZeta_nonneg_all hball) s p

#print axioms smul_factor_eq_split
