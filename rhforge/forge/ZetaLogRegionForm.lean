import LogNormZetaExplicit
import RegionHypForm

/-!
# §1302 — ζ log-growth bound stated in the standard zero-free-region form

Capstone of the §1286–1301 analytic arc: combining the explicit `log‖ζ‖` bound §1300 with the region-shape
lemma §1301, the whole result reads in the classical zero-free-region form — for `s = σ+it` with `0<σ<1`,
`|t| ≥ 1`, and `σ ≥ 1 − 1/log(|t|+1)`:

  `log‖ζ(σ+it)‖ ≤ log(e·log(|t|+1) + 2 + e + e/σ)`.

This is precisely the Borel–Carathéodory ceiling `Re(Log ζ) ≤ M`, `M ≍ log log|t|`, over the region whose
boundary `σ = 1 − 1/log(|t|+1)` matches the shape `dlvp_region_explicit` produces (`β ≤ 1 − 1/(14AL)`).  The
analytic engine of the sharp zero-free region, fully machine-checked and in canonical form.  EXIT0 / std-axioms.
-/

theorem log_norm_riemannZeta_le_in_region {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1)
    (hst : 1 ≤ |s.im|) (hσ : 1 - 1 / Real.log (|s.im| + 1) ≤ s.re) :
    Real.log ‖riemannZeta s‖
      ≤ Real.log (Real.exp 1 * Real.log (|s.im| + 1) + 2 + Real.exp 1 + Real.exp 1 / s.re) :=
  log_norm_riemannZeta_le_log_im hs0 hs1 hst (region_hyp_of_sigma_ge hst hσ)

#print axioms log_norm_riemannZeta_le_in_region
