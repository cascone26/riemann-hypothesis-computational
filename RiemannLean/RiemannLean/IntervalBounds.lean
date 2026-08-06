import Mathlib.Analysis.SpecialFunctions.Exponential

/-!
# Verified interval bounds for the Weil-form entries (the interval-lifting layer)

The one gap in the machine-checked Weil-positivity certificate (LEDGER §755/§764) is a
Lean-verified rigorous bound on the transcendental Weil entries — this would discharge the
`(Q − M).PosSemidef` hypothesis of `weil_posSemidef_of_dominates`.  The Weil entries are
`arch + pole − prime`; the prime and pole terms are *finite* combinations of elementary
functions (`Real.exp`, `Real.cos`) with rational coefficients (`Λ(n)/√n`), and the Gaussian
factor is `e^{−y²/4}`.  This module begins the layer of rigorous rational bounds on those
elementary factors, built on mathlib's `Real.add_one_le_exp` etc.

This is a genuine first brick of a large undertaking (a full Lean-verified interval-arithmetic
evaluator for the transcendental entries, incl. the digamma archimedean integral, does not yet
exist in mathlib — an honest infrastructure gap).
-/

/-- Rigorous rational bound on the Gaussian factor `e^{−y²/4}` at `y² = 1`: `exp(−1/4) ≤ 4/5`.
From `5/4 ≤ exp(1/4)` (`add_one_le_exp`) by taking inverses. -/
theorem exp_neg_quarter_le : Real.exp (-(1 / 4)) ≤ 4 / 5 := by
  have h : (5 : ℝ) / 4 ≤ Real.exp (1 / 4) := by
    have := Real.add_one_le_exp (1 / 4 : ℝ); linarith
  rw [Real.exp_neg, inv_le_comm₀ (Real.exp_pos _) (by norm_num)]
  linarith
