import ZpowToNat
import DivisorLFunctionNonneg

/-!
# §1487 — `LFunction χ`'s `g_full`: ℤ-exponent form = ℕ-exponent form (reconciliation)

L-analogue of ζ §1394.  The half-disk decomposition §1486 carries `g_full` with **ℤ** exponents
`(s-ρ)^{divisor ρ}`, while the factor-bound §1482 uses **ℕ** exponents `(s-ρ)^{(divisor ρ).toNat}`.  Since
`LFunction χ` is entire, `divisor ≥ 0` everywhere (§1479), so `zpow_eq_pow_toNat` makes the two products equal
factor-by-factor — the two `g_full` forms are the *same function*, so their `logDeriv`s agree and §1482's
`Re logDeriv ≥ -B` discharges the `hgB` about §1486's decomposition factor.  Cleaner than ζ (no `hball`).
EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric

theorem gfull_zpow_eq_nat_LFunction {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {T : Finset ℂ} {c : ℂ} {R : ℝ} {g : ℂ → ℂ} :
    (fun s => (∏ ρ ∈ T, (s - ρ) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R) ρ)) * g s)
      = (fun s => (∏ ρ ∈ T,
          (s - ρ) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R) ρ).toNat) * g s) := by
  funext s
  congr 1
  exact Finset.prod_congr rfl
    (fun ρ _ => zpow_eq_pow_toNat (divisor_LFunction_nonneg_all hχ ρ))

#print axioms gfull_zpow_eq_nat_LFunction
