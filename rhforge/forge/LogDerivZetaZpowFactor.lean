import LogDerivZpowMulSplit

/-!
# §1333 — Multiplicity decomposition of `logDeriv ζ` from `ζ =ᶠ P·h` (P the divisor product)

The multiplicity version of §1327: if `ζ` agrees near `z` with `P·h`, `P(s) = ∏_{ρ∈S}(s-ρ)^{m ρ}` (the divisor
product with integer multiplicities `m`) and `h` differentiable nonvanishing at `z`, then
`logDeriv ζ z = ∑_{ρ∈S} m ρ/(z-ρ) + logDeriv h z`.  (logDeriv is local: `logDeriv_congr_nhds` + §1332.)  This is the
exact decomposition the extract_zeros_poles factorization of ζ produces (`P` = divisor product, `h` = the
nonvanishing analytic factor `g`), for zeros of ANY multiplicity — the general `hdecomp` for the region.
EXIT0 / std-axioms.
-/

open Complex Filter Topology

theorem logDeriv_eq_zpow_sum_of_factor {ζ h : ℂ → ℂ} {z : ℂ} {S : Finset ℂ} {m : ℂ → ℤ}
    (hz : ∀ ρ ∈ S, z - ρ ≠ 0) (hhz : h z ≠ 0) (hdh : DifferentiableAt ℂ h z)
    (hfac : ζ =ᶠ[𝓝 z] fun s => (∏ ρ ∈ S, (s - ρ) ^ (m ρ)) * h s) :
    logDeriv ζ z = (∑ ρ ∈ S, (m ρ : ℂ) / (z - ρ)) + logDeriv h z := by
  rw [(logDeriv_congr_nhds hfac).eq_of_nhds, logDeriv_prod_zpow_mul_split hz hhz hdh]

#print axioms logDeriv_eq_zpow_sum_of_factor
