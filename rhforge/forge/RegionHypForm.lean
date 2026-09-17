import ZetaLogImExplicit

/-!
# §1301 — Region hypothesis in standard form: `σ ≥ 1 − 1/log(|t|+1)` ⟹ `(1-σ)·log(|t|+1) ≤ 1`

The §1299/§1300 regime hypothesis `(1-σ)·log(|t|+1) ≤ 1` is exactly the classical zero-free-region shape
`σ ≥ 1 − 1/log(|t|+1)`.  This lemma converts the natural inequality into the product form the bounds consume,
so §1299's log-growth bound reads directly against the region `dlvp_region_explicit` outputs
(`β ≤ 1 − 1/(14AL)`).  Valid for `|t| ≥ 1` (so `log(|t|+1) > 0`).  EXIT0 / std-axioms.
-/

theorem region_hyp_of_sigma_ge {s : ℂ} (hst : 1 ≤ |s.im|)
    (hσ : 1 - 1 / Real.log (|s.im| + 1) ≤ s.re) :
    (1 - s.re) * Real.log (|s.im| + 1) ≤ 1 := by
  have hL : 0 < Real.log (|s.im| + 1) := Real.log_pos (by linarith)
  have hle : 1 - s.re ≤ 1 / Real.log (|s.im| + 1) := by linarith
  calc (1 - s.re) * Real.log (|s.im| + 1)
      ≤ (1 / Real.log (|s.im| + 1)) * Real.log (|s.im| + 1) :=
        mul_le_mul_of_nonneg_right hle hL.le
    _ = 1 := by field_simp

#print axioms region_hyp_of_sigma_ge
