import STOffPairedSum

open scoped MatrixOrder Matrix

theorem S_T_off_sum_eq_ofReal_re (N : ℕ) (L : ℝ) (K : ℕ) (z : Fin K → ℂ)
    (i j : Fin (2 * N + 1)) :
    S_T_off_sum N L K z i j = ((S_T_off_sum N L K z i j).re : ℂ) := by
  have hconj := S_T_off_sum_conj_self N L K z
  have h := congrFun (congrFun hconj i) j
  simp only [Matrix.map_apply] at h
  exact (Complex.conj_eq_iff_re.mp h).symm

#print axioms S_T_off_sum_eq_ofReal_re
