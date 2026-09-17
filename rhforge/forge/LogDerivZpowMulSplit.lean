import LogDerivProdZpow

/-!
# §1332 — Multiplicity product split: `logDeriv((∏(s-ρ)^{m ρ})·h) = ∑ m ρ/(z-ρ) + logDeriv h`

The multiplicity-aware version of §1316 (which used a simple product).  With `P(s) = ∏_{ρ∈S}(s-ρ)^{m ρ}` (the
divisor product, `m : ℂ → ℤ`) and `h` differentiable nonvanishing at `z`, `logDeriv(P·h) = logDeriv P + logDeriv h`
(`logDeriv_mul`) `= ∑ m ρ/(z-ρ) + logDeriv h` (§1329).  This is the shape the extract_zeros_poles factorization
`ζ = P·g` (P the divisor product with integer multiplicities, g nonvanishing) produces — closing the last gap in
the ζ'/ζ decomposition for zeros of any multiplicity.  EXIT0 / std-axioms.
-/

open Complex Finset

theorem logDeriv_prod_zpow_mul_split {S : Finset ℂ} {m : ℂ → ℤ} {h : ℂ → ℂ} {z : ℂ}
    (hz : ∀ ρ ∈ S, z - ρ ≠ 0) (hhz : h z ≠ 0) (hdh : DifferentiableAt ℂ h z) :
    logDeriv (fun s => (∏ ρ ∈ S, (s - ρ) ^ (m ρ)) * h s) z
      = (∑ ρ ∈ S, (m ρ : ℂ) / (z - ρ)) + logDeriv h z := by
  have hPz : (∏ ρ ∈ S, (z - ρ) ^ (m ρ)) ≠ 0 :=
    Finset.prod_ne_zero_iff.mpr (fun ρ hρ => zpow_ne_zero (m ρ) (hz ρ hρ))
  have hdP : DifferentiableAt ℂ (fun s => ∏ ρ ∈ S, (s - ρ) ^ (m ρ)) z :=
    DifferentiableAt.fun_finsetProd (fun ρ hρ => ((differentiableAt_fun_id).sub_const ρ).zpow
      (Or.inl (hz ρ hρ)))
  rw [logDeriv_mul z hPz hhz hdP hdh, logDeriv_prod_zpow_sub hz]

#print axioms logDeriv_prod_zpow_mul_split
