import STOffRealSymm
import ContrastMatrixPacked

open scoped MatrixOrder Matrix

/-!
# Off-critical-line lane: `G^{off,0}_{N,L} := Cᵀ·S^{(T,off)}·C` is symmetric (Shi Eq 7.14)

Assembles the paper's actual Eq (7.14) compressed matrix `G^{off,0}_{N,L}` from §2359's real,
symmetric `S_T_off_real` and §2335's already-built, generic, injective-`mulVec` `C_matrix` — the
first genuine connection between Section 7's concrete off-critical-line construction and the
Section-4-lane machinery built much earlier this session. Note: unlike Section 4/9's `G_{N,L}`
(built from the always-PosSemidef shifted Weil matrix `W_NL`), Section 7's `S^{(T,off)}` is
"generally indefinite" (Eq 7.3, the paper's own words: "It has no fixed sign") — so `G^{off,0}`
CANNOT be shown PosSemidef, only symmetric (a standard congruence fact independent of sign:
`(Cᵀ S C)ᵀ = Cᵀ Sᵀ C = Cᵀ S C` whenever `S` is symmetric). Search-before-build: grepped for
`G_off0`/`GOff0`/`K_off0`/`KOff0` — fresh. Fully general, no zeta/RH dependence. -/
theorem congruence_isSymm {n m : Type*} [Fintype n] [Fintype m]
    (S : Matrix n n ℝ) (hS : S.IsSymm) (C : Matrix n m ℝ) :
    (Cᵀ * S * C).IsSymm := by
  unfold Matrix.IsSymm at hS ⊢
  rw [Matrix.transpose_mul, Matrix.transpose_mul, Matrix.transpose_transpose, hS, Matrix.mul_assoc]

noncomputable def G_off0 (N : ℕ) (L : ℝ) (K : ℕ) (z : Fin K → ℂ) :
    Matrix (Fin (2 * N)) (Fin (2 * N)) ℝ :=
  (C_matrix (2 * N + 1) (by omega))ᵀ * (S_T_off_real N L K z) * (C_matrix (2 * N + 1) (by omega))

theorem G_off0_isSymm (N : ℕ) (L : ℝ) (K : ℕ) (z : Fin K → ℂ) :
    (G_off0 N L K z).IsSymm :=
  congruence_isSymm (S_T_off_real N L K z) (S_T_off_real_isSymm N L K z)
    (C_matrix (2 * N + 1) (by omega))

#print axioms congruence_isSymm
#print axioms G_off0_isSymm
