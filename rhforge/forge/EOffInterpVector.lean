import P2KLDerivAtNode
import QOffKPoly

open scoped MatrixOrder Matrix

/-!
# Off-critical-line lane: Eq (7.13)'s interpolation coefficients `c^{off}`/`e^{off}`

Shi Eq (7.13): `c^{off}_{n,2K,L} := Q^{off}_K(ν_{n,L})/P'_{2K,L}(ν_{n,L})`, `e^{off}_{n,2K,L} :=
√L·c^{off}_{n,2K,L}`. This assembles EVERY piece built across §2364-2368 into the paper's actual
explicit interpolation coefficient: `Q_off_K` at each frequency node, divided by `P_2KL`'s
derivative there (using §2368's `P_2KL_hasDerivAt_node`, converted to `deriv`-notation via
`HasDerivAt.deriv`, matching the paper's own `P'` notation literally). Search-before-build: grepped
for `c_off`/`e_off`/`COff`/`EOff` — fresh (two coincidental substring hits inside unrelated
`_offline` identifiers, confirmed harmless). Fully general, no zeta/RH dependence. -/
noncomputable def c_off (K : ℕ) (L : ℝ) (z : Fin K → ℂ) (n : Fin (2 * (2 * K) + 1)) : ℝ :=
  Q_off_K K z (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L)
    / (deriv (P_2KL K L) (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L))

noncomputable def e_off (K : ℕ) (L : ℝ) (z : Fin K → ℂ) (n : Fin (2 * (2 * K) + 1)) : ℝ :=
  Real.sqrt L * c_off K L z n

theorem c_off_eq_div_prod (K : ℕ) (L : ℝ) (z : Fin K → ℂ) (n : Fin (2 * (2 * K) + 1)) :
    c_off K L z n = Q_off_K K z (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L)
      / (∏ j ∈ Finset.univ.erase n,
          ((2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L)
            - (2 * Real.pi * (((j : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L))) := by
  unfold c_off
  rw [(P_2KL_hasDerivAt_node K L n).deriv]

#print axioms c_off_eq_div_prod
