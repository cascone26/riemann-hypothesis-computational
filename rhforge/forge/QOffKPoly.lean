import Mathlib

open scoped MatrixOrder Matrix

/-!
# Off-critical-line lane: `Q^{off}_K(w)` (Eq 7.11-7.12), the quartet interpolation polynomial

Shi Eq (7.11): `Q^{off}_K(w) := ∏_{k=1}^K (w²−z_k²)(w²−z̄_k²)`, with the stated real-coefficient
expansion (Eq 7.12) `(w²−z_k²)(w²−z̄_k²) = w⁴ − 2(γ_k²−η_k²)w² + (γ_k²+η_k²)²` for `z_k=γ_k+iη_k`.
This brick defines BOTH forms — the abstract complex product (`Q_off_K_complex`, matching Eq 7.11
literally) and the concrete real quartic-product closed form (`Q_off_K`, matching Eq 7.12) — and
proves they agree, i.e. the complex product IS genuinely real and equals this explicit real
polynomial value. This is the actual interpolation-polynomial machinery Theorem 7.2's proof needs
(reusing, per the paper's own statement, "the partial-fraction argument of Theorem 5.2" — whose
general abstract version already exists in the corpus as `rational_function_partial_fraction`,
confirmed via search-before-build). Search-before-build: grepped for `P_2K_L`/`Q_off_K`/`QoffK`/
`P2KL` — fresh. Fully general, no zeta/RH dependence beyond the standing complex-parameter setup.
-/
noncomputable def Q_off_K_complex (K : ℕ) (z : Fin K → ℂ) (w : ℂ) : ℂ :=
  ∏ k : Fin K, ((w ^ 2 - (z k) ^ 2) * (w ^ 2 - (starRingEnd ℂ (z k)) ^ 2))

noncomputable def Q_off_K (K : ℕ) (z : Fin K → ℂ) (w : ℝ) : ℝ :=
  ∏ k : Fin K, (w ^ 4 - 2 * ((z k).re ^ 2 - (z k).im ^ 2) * w ^ 2
    + ((z k).re ^ 2 + (z k).im ^ 2) ^ 2)

theorem Q_off_K_complex_eq_ofReal (K : ℕ) (z : Fin K → ℂ) (w : ℝ) :
    Q_off_K_complex K z (w : ℂ) = ((Q_off_K K z w : ℝ) : ℂ) := by
  unfold Q_off_K_complex Q_off_K
  push_cast
  apply Finset.prod_congr rfl
  intro k _
  have hzk : z k = Complex.mk (z k).re (z k).im := Complex.ext rfl rfl
  apply Complex.ext
  · rw [hzk]
    simp [Complex.mul_re, Complex.sub_re, Complex.sub_im, pow_two, Complex.mul_im,
      Complex.conj_re, Complex.conj_im, ← Complex.ofReal_pow]
    ring
  · rw [hzk]
    simp [Complex.mul_im, Complex.sub_re, Complex.sub_im, pow_two, Complex.mul_re,
      Complex.conj_re, Complex.conj_im, ← Complex.ofReal_pow]
    ring

theorem Q_off_K_even (K : ℕ) (z : Fin K → ℂ) (w : ℝ) :
    Q_off_K K z (-w) = Q_off_K K z w := by
  unfold Q_off_K
  apply Finset.prod_congr rfl
  intro k _
  ring

#print axioms Q_off_K_complex_eq_ofReal
#print axioms Q_off_K_even
