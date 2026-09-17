import Mathlib.Analysis.Calculus.LogDeriv
import Mathlib.Tactic

set_option maxHeartbeats 1000000

/-!
# §1676 — `logDeriv` of the zeros polynomial with multiplicities: `Σ m/(z-ρ)`

The ℕ-power version of the tower's `logDeriv_prod_sub` (§ LogDerivProdSub, multiplicities 1):
`logDeriv (∏_ρ (s-ρ)^{m ρ}) z = Σ_ρ m ρ/(z-ρ)` away from the zeros. The zeros' partial-fraction
term of the everywhere-factorized Landau split (§1677).
-/

theorem logDeriv_pow_sub {ρ z : ℂ} (m : ℕ) (hz : z - ρ ≠ 0) :
    logDeriv (fun s : ℂ => (s - ρ) ^ m) z = (m : ℂ) / (z - ρ) := by
  cases m with
  | zero =>
    simp [logDeriv_const]
  | succ n =>
    have hd : HasDerivAt (fun s : ℂ => (s - ρ) ^ (n + 1))
        ((n + 1 : ℕ) * (z - ρ) ^ n * 1) z :=
      ((hasDerivAt_id z).sub_const ρ).pow (n + 1)
    rw [logDeriv_apply, hd.deriv]
    have hpne : (z - ρ) ^ n ≠ 0 := pow_ne_zero n hz
    rw [pow_succ]
    field_simp

theorem logDeriv_prod_pow {S : Finset ℂ} {m : ℂ → ℕ} {z : ℂ} (hz : ∀ ρ ∈ S, z - ρ ≠ 0) :
    logDeriv (fun s : ℂ => ∏ ρ ∈ S, (s - ρ) ^ (m ρ)) z = ∑ ρ ∈ S, (m ρ : ℂ) / (z - ρ) := by
  have h := logDeriv_prod (s := S) (f := fun ρ => fun s : ℂ => (s - ρ) ^ (m ρ)) (x := z)
    (fun ρ hρ => pow_ne_zero _ (hz ρ hρ))
    (fun ρ _ => (differentiableAt_id.sub_const ρ).pow _)
  rw [show (fun s : ℂ => ∏ ρ ∈ S, (s - ρ) ^ (m ρ))
      = (∏ ρ ∈ S, (fun s : ℂ => (s - ρ) ^ (m ρ)) ·) from ?_, h]
  · exact Finset.sum_congr rfl (fun ρ hρ => logDeriv_pow_sub (m ρ) (hz ρ hρ))
  · funext s
    simp

#print axioms logDeriv_pow_sub
#print axioms logDeriv_prod_pow
