import LogDerivZetaDecomp

/-!
# §1327 — Decomposition of `logDeriv ζ` from a local factorization `ζ = P·h`

If `ζ` agrees near `z` with `P·h`, `P(s) = ∏_{ρ∈S}(s-ρ)` (the disk zeros) and `h` differentiable nonvanishing at
`z`, then `logDeriv ζ z = ∑_{ρ∈S} 1/(z-ρ) + logDeriv h z`.  (Logarithmic derivative is local — `logDeriv_congr_nhds`
— so it only sees the factorization on a neighbourhood.)  This reduces the master-composition hypothesis `hdecomp`
(§1326) to the single fact `ζ =ᶠ P·h`, i.e. that `P` carries exactly ζ's disk zeros (with multiplicity) and
`h = ζ/P` is the analytic nonvanishing quotient.  EXIT0 / std-axioms.
-/

open Complex Filter Topology

theorem logDeriv_eq_sum_add_of_factor {ζ h : ℂ → ℂ} {z : ℂ} {S : Finset ℂ}
    (hz : ∀ ρ ∈ S, z - ρ ≠ 0) (hhz : h z ≠ 0) (hdh : DifferentiableAt ℂ h z)
    (hfac : ζ =ᶠ[𝓝 z] fun s => (∏ ρ ∈ S, (s - ρ)) * h s) :
    logDeriv ζ z = (∑ ρ ∈ S, 1 / (z - ρ)) + logDeriv h z := by
  rw [(logDeriv_congr_nhds hfac).eq_of_nhds, logDeriv_prod_mul_split hz hhz hdh]

#print axioms logDeriv_eq_sum_add_of_factor
