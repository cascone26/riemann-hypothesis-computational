import LogDerivProdSub

/-!
# §1316 — Product split of the log-derivative: `logDeriv(P·h) = ∑1/(s-ρ) + logDeriv h`

Writing `ζ = P·h` on a disk, `P(s) = ∏_{ρ∈S}(s-ρ)` (the finite disk-zeros) and `h = ζ/P` the nonvanishing
quotient, the log-derivative splits (Mathlib `logDeriv_mul`): `logDeriv ζ = logDeriv P + logDeriv h`.  With §1315
(`logDeriv P = ∑1/(s-ρ)`) this is EXACTLY the decomposition `F = (∑1/(s-ρ)) + g` that §1314 turns into `hzero`,
with `g = logDeriv h`.  So the ENTIRE quantitative zero-free region is now algebraically reduced to a single
analytic estimate: `Re(logDeriv h) ≥ -B` on the disk (the Borel–Carathéodory bound on the nonvanishing quotient).
EXIT0 / std-axioms.
-/

open Complex

theorem logDeriv_prod_mul_split {S : Finset ℂ} {h : ℂ → ℂ} {z : ℂ}
    (hz : ∀ ρ ∈ S, z - ρ ≠ 0) (hhz : h z ≠ 0) (hdh : DifferentiableAt ℂ h z) :
    logDeriv (fun s => (∏ ρ ∈ S, (s - ρ)) * h s) z = (∑ ρ ∈ S, 1 / (z - ρ)) + logDeriv h z := by
  have hPz : (∏ ρ ∈ S, (z - ρ)) ≠ 0 := Finset.prod_ne_zero_iff.mpr hz
  have hdP : DifferentiableAt ℂ (fun s => ∏ ρ ∈ S, (s - ρ)) z :=
    DifferentiableAt.fun_finsetProd (fun ρ _ => (differentiableAt_fun_id).sub_const ρ)
  rw [logDeriv_mul z hPz hhz hdP hdh, logDeriv_prod_sub hz]

#print axioms logDeriv_prod_mul_split
