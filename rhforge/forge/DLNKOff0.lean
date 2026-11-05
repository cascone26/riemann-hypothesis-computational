import GOff0Symm

open scoped MatrixOrder Matrix

/-!
# Off-critical-line lane: `D_{L,N}` (frequency diagonal matrix) + `K^{off,0}_{N,L}` (Eq 7.15)

Search-before-build: grepped for `D_LN`/`D_L_N`/`DLN`/`freqMatrix`/`diagonal.*nu_n` — genuinely
fresh. `D_{L,N} := diag(ν_{n,L})_{n∈I_N}` is used throughout Shi's paper (Section 4's `K_{N,L} :=
C*WD_{L,N}C`, Eq 4.9, and Section 7's `K^{off,0}` here) but had never actually been built anywhere
in this session's Shi-lane corpus — the Section 4 lane's own `K_{N,L}` construction was apparently
deferred at the time and never revisited before this session moved on to other threads. `D_LN`
closes that standing gap generally (usable by BOTH Section 4's `K_{N,L}` and Section 7's
`K^{off,0}` alike), matching the established `ν_{n,L}=2π(n-N)/L` indexing convention used
throughout this session (`S_NL`, `vPlusL`, etc.). `K_off0` then assembles Shi's actual Eq (7.15)
object. Fully general, no zeta/RH dependence. -/
noncomputable def D_LN (N : ℕ) (L : ℝ) : Matrix (Fin (2 * N + 1)) (Fin (2 * N + 1)) ℝ :=
  Matrix.diagonal (fun n : Fin (2 * N + 1) => 2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L)

theorem D_LN_isSymm (N : ℕ) (L : ℝ) : (D_LN N L).IsSymm := by
  unfold D_LN Matrix.IsSymm
  exact Matrix.diagonal_transpose _

noncomputable def K_off0 (N : ℕ) (L : ℝ) (K : ℕ) (z : Fin K → ℂ) :
    Matrix (Fin (2 * N)) (Fin (2 * N)) ℝ :=
  (C_matrix (2 * N + 1) (by omega))ᵀ * (S_T_off_real N L K z * D_LN N L)
    * (C_matrix (2 * N + 1) (by omega))

#print axioms D_LN_isSymm
