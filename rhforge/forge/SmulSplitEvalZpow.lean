import SmulFactorizedEval

/-!
# §1392 — ℤ-exponent split-evaluation (for the half-disk logDeriv decomposition)

The `ℤ`-zpow analogue of §1386 (no `toNat`), matching §1336's `m : ℂ → ℤ = divisor` and §1332's zpow product.
`((∏ᶠ (·-u)^{d u}) • g) s = P_half(s) · (P_rest(s) · g s)` with `ℤ` exponents, via §1338 (smul-eval) and the
Finset filter split.  This feeds the **half-disk decomposition** `logDeriv ζ z = ∑_{half-disk} divisor/(z-ρ) +
logDeriv g_full z` (the corrected §1349-analogue: keep only the half-disk zeros in the sum, the annulus zeros live
in `g_full` and are bounded by the two-radii §1378).  EXIT0 / std-axioms.
-/

open Function

open Classical in
theorem smul_split_eval_zpow {d : ℂ → ℤ} (hd : (support d).Finite) {g : ℂ → ℂ} (s : ℂ)
    (p : ℂ → Prop) :
    ((∏ᶠ u, (· - u) ^ (d u)) • g) s
      = (∏ u ∈ hd.toFinset.filter p, (s - u) ^ (d u))
        * ((∏ u ∈ hd.toFinset.filter (fun u => ¬ p u), (s - u) ^ (d u)) * g s) := by
  rw [smul_factorized_rational_eval hd, ← Finset.prod_filter_mul_prod_filter_not hd.toFinset p]
  ring

#print axioms smul_split_eval_zpow
