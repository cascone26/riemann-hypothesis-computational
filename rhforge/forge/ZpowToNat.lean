import Mathlib.Analysis.SpecialFunctions.Complex.Circle

/-!
# §1380 — `a^n = a^(n.toNat)` for `0 ≤ n` (zpow→nat bridge for the factorization)

`extract_zeros_poles` produces the product `∏ᶠ u, (·-u)^(divisor u)` with **ℤ** exponents (`zpow`), but the
factor-bound capstone §1378 consumes a product with **ℕ** exponents (`m : ℂ → ℕ`).  On a pole-free ball the
divisor is `≥ 0` (§1365), so each `(·-u)^(divisor u : ℤ) = (·-u)^((divisor u).toNat)`.  This is the elementary
bridge `a^n = a^(n.toNat)` for `0 ≤ n`, letting `m = fun u => (divisor u).toNat` line the two products up.
EXIT0 / std-axioms.
-/

theorem zpow_eq_pow_toNat {a : ℂ} {n : ℤ} (hn : 0 ≤ n) : a ^ n = a ^ n.toNat := by
  conv_lhs => rw [← Int.toNat_of_nonneg hn]
  rw [zpow_natCast]

#print axioms zpow_eq_pow_toNat
