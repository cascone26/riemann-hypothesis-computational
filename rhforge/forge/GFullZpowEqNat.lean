import ZpowToNat
import DivisorNonnegAll

/-!
# §1394 — `g_full` ℤ-exponent form = ℕ-exponent form (reconciliation)

§1393's half-disk decomposition carries `g_full` with **ℤ** exponents `(s-ρ)^{divisor ρ}` (matching §1336/§1332),
while §1391/§1378's factor-bound uses **ℕ** exponents `(s-ρ)^{(divisor ρ).toNat}` (matching the norm chain).  On a
pole-free ball `divisor ≥ 0` (§1387), so `zpow_eq_pow_toNat` (§1380) makes the two products equal factor-by-factor
— hence the two `g_full` forms are the *same function*, so their `logDeriv`s agree and §1378's `Re logDeriv ≥ -B`
discharges §1336's `hgB` about §1393's decomposition factor.  EXIT0 / std-axioms.
-/

open Metric

theorem gfull_zpow_eq_nat {T : Finset ℂ} {c : ℂ} {R : ℝ} (hball : (1 : ℂ) ∉ ball c R) {g : ℂ → ℂ} :
    (fun s => (∏ ρ ∈ T, (s - ρ) ^ (MeromorphicOn.divisor riemannZeta (ball c R) ρ)) * g s)
      = (fun s => (∏ ρ ∈ T, (s - ρ) ^ (MeromorphicOn.divisor riemannZeta (ball c R) ρ).toNat) * g s) := by
  funext s
  congr 1
  exact Finset.prod_congr rfl
    (fun ρ _ => zpow_eq_pow_toNat (divisor_riemannZeta_nonneg_all hball ρ))

#print axioms gfull_zpow_eq_nat
