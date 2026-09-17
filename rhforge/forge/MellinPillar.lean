import TruncatedMellin
import BoundaryTendstoZero
import SumIocTendstoLSeries
import PsiCpowIntegrableOnIoi
import LogDerivZetaVonMangoldt
import Mathlib.MeasureTheory.Integral.IntegralEqImproper

set_option maxHeartbeats 1000000

/-!
# §1588 ★★★★★ — THE MELLIN PILLAR: `-ζ'/ζ(s) = s·∫₁^∞ ψ(t)·t^{-s-1} dt` for `Re s > 1`

The first pillar of the explicit formula.  Take `N → ∞` in the truncated representation §1584: the boundary term
dies (§1585), the truncated sums converge to `L(Λ) = -ζ'/ζ` (§1586, §1356), and the truncated integrals converge
to the improper integral (`intervalIntegral_tendsto_integral_Ioi` with §1587).  Limit uniqueness closes it.
This identity is what contour integration will turn into the explicit formula
`ψ(x) = x - ∑_ρ x^ρ/ρ - …`, where our zero-free regions (§1575/§1577) bound the zero sum — the road to
`ψ(x) = x + O(x·exp(-c√log x))`.  EXIT0 / std-axioms.
-/

open Complex Chebyshev ArithmeticFunction Filter MeasureTheory
open scoped LSeries.notation

theorem neg_logDeriv_zeta_eq_mellin {s : ℂ} (hs : 1 < s.re) :
    -logDeriv riemannZeta s
      = s * ∫ t : ℝ in Set.Ioi (1 : ℝ), ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ (-s - 1) := by
  have hs0 : s ≠ 0 := by
    intro h; rw [h] at hs; norm_num [Complex.zero_re] at hs
  have hσ0 : (0 : ℝ) ≤ s.re := by linarith
  -- LHS limit
  have hL := sum_Ioc_tendsto_LSeries hs
  -- RHS limit
  have hI := MeasureTheory.intervalIntegral_tendsto_integral_Ioi (1 : ℝ)
    (psi_cpow_integrableOn_Ioi hs) tendsto_natCast_atTop_atTop
  have hB := boundary_term_tendsto_zero hs
  have hR : Tendsto (fun N : ℕ => (N : ℂ) ^ (-s) * ((ψ (N : ℝ) : ℝ) : ℂ)
      + s * ∫ t : ℝ in (1 : ℝ)..(N : ℝ), ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ (-s - 1))
      atTop (nhds (0 + s * ∫ t : ℝ in Set.Ioi (1 : ℝ), ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ (-s - 1))) :=
    hB.add (hI.const_mul s)
  -- the two sequences agree eventually (N ≥ 1)
  have hcongr : Tendsto (fun N : ℕ => ∑ n ∈ Finset.Ioc 0 N, (Λ n : ℂ) * (n : ℂ) ^ (-s))
      atTop (nhds (0 + s * ∫ t : ℝ in Set.Ioi (1 : ℝ), ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ (-s - 1))) := by
    apply hR.congr'
    filter_upwards [eventually_ge_atTop 1] with N hN
    exact (truncated_mellin hs0 hσ0 hN).symm
  have huniq := tendsto_nhds_unique hL hcongr
  rw [logDeriv_riemannZeta_eq_neg_vonMangoldt hs, neg_neg, huniq, zero_add]

#print axioms neg_logDeriv_zeta_eq_mellin
