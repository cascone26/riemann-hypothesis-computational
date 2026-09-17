import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.NumberTheory.Divisors
import Mathlib.Data.Nat.PrimeFin
import Mathlib.Data.Nat.Factorization.Basic

set_option maxHeartbeats 400000

/-!
# §1562 — `∑_{p|N} log p ≤ log N`  (the finite-product correction is at most `log N`)

Helper for the QUADRATIC capstones: the correction constant `∑_{p∈N.primeFactors} log p` (from §1554/§1556) is
bounded by `log N`, since `∏_{p|N} p ∣ N` (`Nat.prod_primeFactors_dvd`) so `∏_{p|N} p ≤ N`, and `∑ log p =
log(∏ p)` (`Real.log_prod`).  This lets the fixed correction be absorbed into the region threshold for large `|γ|`.
EXIT0 / std-axioms.
-/

theorem sum_log_primeFactors_le_log {N : ℕ} (hN : 1 ≤ N) :
    ∑ p ∈ N.primeFactors, Real.log p ≤ Real.log N := by
  have hne : ∀ p ∈ N.primeFactors, (p : ℝ) ≠ 0 := by
    intro p hp
    have := (Nat.prime_of_mem_primeFactors hp).two_le
    exact_mod_cast (by omega : p ≠ 0)
  have heq : ∑ p ∈ N.primeFactors, Real.log (p : ℝ)
      = Real.log (∏ p ∈ N.primeFactors, (p : ℝ)) :=
    (Real.log_prod hne).symm
  rw [heq]
  have hdvd : (∏ p ∈ N.primeFactors, p) ∣ N := Nat.prod_primeFactors_dvd N
  have hle : (∏ p ∈ N.primeFactors, p) ≤ N := Nat.le_of_dvd (by omega) hdvd
  have hpos : 0 < ∏ p ∈ N.primeFactors, p :=
    Finset.prod_pos (fun p hp => (Nat.prime_of_mem_primeFactors hp).pos)
  have hcast : (∏ p ∈ N.primeFactors, (p : ℝ)) = ((∏ p ∈ N.primeFactors, p : ℕ) : ℝ) := by
    push_cast; ring
  rw [hcast]
  apply Real.log_le_log (by exact_mod_cast hpos)
  exact_mod_cast hle

#print axioms sum_log_primeFactors_le_log
