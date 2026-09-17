import FinprodZpowFinset
import Mathlib.Analysis.Meromorphic.FactorizedRational

/-!
# §1338 — Pointwise value of the extract_zeros_poles factor: `((∏ᶠ(·-u)^{d u})•g) s = (∏_{support}(s-u)^{d u})·g s`

`(φ•g) s = φ s · g s` (`Pi.smul_apply'`+`smul_eq_mul`); the finprod-of-functions evaluates pointwise `φ s = ∏ᶠ u,(s-u)^{d u}`
(`finprod_apply`, whose finite-mul-support hypothesis follows from `support d` finite since `d u = 0 ⟹` the factor is `1`);
then §1337 collapses onto the divisor support.  The `smul→mul`+evaluation converter (#2) of the extract_zeros_poles→ζ
integration.  EXIT0 / std-axioms.
-/

open Function

theorem smul_factorized_rational_eval {d : ℂ → ℤ} (hd : (Function.support d).Finite) {g : ℂ → ℂ} (s : ℂ) :
    ((∏ᶠ u, (· - u) ^ (d u)) • g) s = (∏ u ∈ hd.toFinset, (s - u) ^ (d u)) * g s := by
  have hfms : HasFiniteMulSupport (fun u : ℂ => (· - u) ^ (d u)) := by
    apply hd.subset
    intro u hu
    rw [Function.mem_support]
    intro hdu
    exact hu (by ext x; simp [hdu])
  rw [Pi.smul_apply', smul_eq_mul, finprod_apply hfms]
  simp only [Pi.pow_apply]
  rw [finprod_zpow_sub_eq_finset_prod hd]

#print axioms smul_factorized_rational_eval
