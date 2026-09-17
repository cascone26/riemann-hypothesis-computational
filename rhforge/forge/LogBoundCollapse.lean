import Mathlib.Analysis.SpecialFunctions.Log.Basic

/-!
# §1418 — log collapse `log(M/P) - log G ≤ log(M/Z)` when `Z ≤ P·G`

Pure real lemma: the `P`-term (which will be `(R/2)^N`) cancels.  With `Z ≤ P·G` and all positive,
`log Z ≤ log(P·G) = log P + log G`, so `log(M/P) - log G = log M - log P - log G ≤ log M - log Z = log(M/Z)`.
This is the algebra that collapses §1402's `BOUND = 4(log(Mζ/(R/2)^N) - log‖g_full c‖)/(R/2)` to the explicit
`4·log(Mζ/‖ζ c‖)/(R/2)` (`Z = ‖ζ c‖`, `P = (R/2)^N`, `G = ‖g_full c‖`, `M = Mζ`).  EXIT0 / std-axioms.
-/

theorem log_bound_collapse {M P G Z : ℝ} (hM : 0 < M) (hP : 0 < P) (hG : 0 < G) (hZ : 0 < Z)
    (h : Z ≤ P * G) : Real.log (M / P) - Real.log G ≤ Real.log (M / Z) := by
  rw [Real.log_div hM.ne' hP.ne', Real.log_div hM.ne' hZ.ne']
  have hlog : Real.log Z ≤ Real.log P + Real.log G := by
    rw [← Real.log_mul hP.ne' hG.ne']
    exact Real.log_le_log hZ h
  linarith

#print axioms log_bound_collapse
