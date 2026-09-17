import ZetaLogImExplicit

/-!
# §1300 — Explicit `log‖ζ(σ+it)‖ ≤ log(e·log(|t|+1)+2+e+e/σ)` (the Borel–Carathéodory `M` input)

Taking `Real.log` of the |t|-explicit bound §1299.  `Re (Log ζ(s)) = log‖ζ(s)‖`, so this is a pointwise explicit
ceiling `M(t,σ)` on the real part of `log ζ` — exactly what `deriv_bound_of_re_le_general` (§1237) consumes to
turn a bounded real part into a bounded log-derivative `‖ζ'/ζ‖`.  It replaces the forge's *qualitative*
`exists_norm_bound_riemannZeta_sphere` (existential `C`) with a *quantitative* `M ≍ log log|t|`.  Valid for
`0<σ<1`, `|t|≥1`, region `(1-σ)·log(|t|+1) ≤ 1`.  EXIT0 / std-axioms.
-/

theorem log_norm_riemannZeta_le_log_im {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1)
    (hst : 1 ≤ |s.im|) (hreg : (1 - s.re) * Real.log (|s.im| + 1) ≤ 1) :
    Real.log ‖riemannZeta s‖
      ≤ Real.log (Real.exp 1 * Real.log (|s.im| + 1) + 2 + Real.exp 1 + Real.exp 1 / s.re) := by
  have hbound := norm_riemannZeta_le_log_im hs0 hs1 hst hreg
  have hlognn : 0 ≤ Real.log (|s.im| + 1) := Real.log_nonneg (by linarith)
  have hemul : 0 ≤ Real.exp 1 * Real.log (|s.im| + 1) := mul_nonneg (Real.exp_pos 1).le hlognn
  have hediv : 0 ≤ Real.exp 1 / s.re := (div_pos (Real.exp_pos 1) hs0).le
  have hexp : 0 ≤ Real.exp 1 := (Real.exp_pos 1).le
  have hRHSge1 : 1 ≤ Real.exp 1 * Real.log (|s.im| + 1) + 2 + Real.exp 1 + Real.exp 1 / s.re := by
    linarith
  rcases eq_or_lt_of_le (norm_nonneg (riemannZeta s)) with h0 | hpos
  · rw [← h0, Real.log_zero]; exact Real.log_nonneg hRHSge1
  · exact Real.log_le_log hpos hbound

#print axioms log_norm_riemannZeta_le_log_im
