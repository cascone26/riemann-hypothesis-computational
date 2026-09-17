import Mathlib.Analysis.Complex.JensenFormula
import Mathlib.NumberTheory.LSeries.RiemannZeta
import ExplicitFormula

/-!
# §1362 — `circleAverage(log‖ζ‖) c R ≤ log M` from a circle bound `‖ζ‖ ≤ M`

The second brick of the Jensen/factor-bound sub-arc.  If `‖ζ‖ ≤ M` (`M ≥ 1`) on the circle `sphere c |R|`,
then `log‖ζ‖ ≤ log M` pointwise there (the junk value `log 0 = 0 ≤ log M` covers the zeros), so by
`circleAverage_mono_on_of_le_circle` the circle average is `≤ log M`.  Combined with §1361 (Jensen) this
turns the ζ-growth bound (§1130/§1299) into an upper bound on the divisor-weighted zero count — the route to
`‖ζ/P‖ ≤ M'` and §1355's `hgB`.  EXIT0 / std-axioms.
-/

open Complex Metric Real

theorem circleAverage_log_norm_riemannZeta_le {c : ℂ} {R M : ℝ} (hM : 1 ≤ M)
    (hbound : ∀ s ∈ sphere c |R|, ‖riemannZeta s‖ ≤ M) :
    circleAverage (fun s => Real.log ‖riemannZeta s‖) c R ≤ Real.log M := by
  have hci : CircleIntegrable (fun s => Real.log ‖riemannZeta s‖) c R :=
    (meromorphicOn_riemannZeta.mono_set (Set.subset_univ _)).circleIntegrable_log_norm
  apply circleAverage_mono_on_of_le_circle hci
  intro s hs
  rcases eq_or_lt_of_le (norm_nonneg (riemannZeta s)) with h0 | hpos
  · rw [← h0, Real.log_zero]; exact Real.log_nonneg hM
  · exact Real.log_le_log hpos (hbound s hs)

#print axioms circleAverage_log_norm_riemannZeta_le
