import STOffRealValued

open scoped MatrixOrder Matrix

/-!
# Off-critical-line lane: the real matrix `S^{(T,off)}_{N,L}` is symmetric (Shi Section 7.1 opening claim)

Defines the actual real-valued `S_T_off_real` (§2354-2355's entrywise real-part bridge, packaged as
a genuine `Matrix _ _ ℝ`), and proves it is symmetric — closing the FIRST half of Section 7.1's own
opening sentence: "`S^{(T,off)}_{N,L}` is real symmetric." (The real-valuedness half was §2354-2355;
this brick closes symmetry.) Route: `H_L_G_C`'s own building blocks are `vecMulVec f f` terms, which
are ALWAYS symmetric (`(vecMulVec a a) i j = a i * a j = a j * a i` via `mul_comm`, no injectivity/
positivity needed) — this symmetry survives entrywise real-part extraction and finite sums
unconditionally. Search-before-build: grepped for `S_T_off_real`/`STOffReal` — fresh. Fully general,
no zeta/RH dependence. -/
theorem vecMulVec_isSymm_self {n : Type*} (a : n → ℂ) : (Matrix.vecMulVec a a).IsSymm := by
  unfold Matrix.IsSymm
  ext i j
  simp [Matrix.vecMulVec, Matrix.transpose_apply, mul_comm]

theorem H_L_G_C_isSymm (N : ℕ) (L : ℝ) (z : ℂ) : (H_L_G_C N L z).IsSymm := by
  unfold H_L_G_C
  exact (vecMulVec_isSymm_self (vPlusLC N L z)).add (vecMulVec_isSymm_self (vMinusLC N L z))

theorem S_T_off_sum_isSymm (N : ℕ) (L : ℝ) (K : ℕ) (z : Fin K → ℂ) :
    (S_T_off_sum N L K z).IsSymm := by
  unfold S_T_off_sum Matrix.IsSymm
  rw [Matrix.transpose_sum]
  apply Finset.sum_congr rfl
  intro k _
  have h1 := (H_L_G_C_isSymm N L (z k))
  have h2 := (H_L_G_C_isSymm N L (starRingEnd ℂ (z k)))
  unfold Matrix.IsSymm at h1 h2
  rw [Matrix.transpose_add, h1, h2]

noncomputable def S_T_off_real (N : ℕ) (L : ℝ) (K : ℕ) (z : Fin K → ℂ) :
    Matrix (Fin (2 * N + 1)) (Fin (2 * N + 1)) ℝ :=
  fun i j => (S_T_off_sum N L K z i j).re

theorem S_T_off_real_isSymm (N : ℕ) (L : ℝ) (K : ℕ) (z : Fin K → ℂ) :
    (S_T_off_real N L K z).IsSymm := by
  unfold S_T_off_real Matrix.IsSymm
  funext i j
  show (S_T_off_sum N L K z j i).re = (S_T_off_sum N L K z i j).re
  have h := congrFun (congrFun (S_T_off_sum_isSymm N L K z) i) j
  simp only [Matrix.transpose_apply] at h
  rw [h]

#print axioms H_L_G_C_isSymm
#print axioms S_T_off_sum_isSymm
#print axioms S_T_off_real_isSymm
