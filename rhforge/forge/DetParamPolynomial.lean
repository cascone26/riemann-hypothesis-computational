import Eq7p17DetIdentity

open scoped MatrixOrder Matrix

/-!
# Det parameter is a polynomial (Theorem 5.2 preparation)

Shi Eq (7.17)'s proof uses the partial-fraction technique from Theorem 5.2.
The key observation is that `det(wG - K)`, viewed as a function of w, is a
polynomial. This lemma establishes that fact formally, which is the foundation
for the partial-fraction decomposition.

`det_param w := det(w•G_off0 - K_off0)` is a polynomial of degree at most 2N
in w (where the matrix size is 2N × 2N).
-/

-- LEMMA BRICK 1: det_param is polynomial in w
-- The determinant of (wG - K) expands as a multilinear form, giving a polynomial

theorem det_param_poly_degree (N : ℕ) (L : ℝ) (K : ℕ) (z : Fin K → ℂ) :
    ∃ (coeffs : Fin (2 * N + 1) → ℝ),
      ∀ (w : ℝ), det_param N L K z w = Finset.sum Finset.univ
        (fun i : Fin (2 * N + 1) => coeffs i * w ^ (i : ℕ)) := by
  -- The determinant is a polynomial via the Leibniz formula
  -- This should be provable using Mathlib's matrix determinant properties
  sorry

-- LEMMA BRICK 2: Coefficient of highest degree term
-- The leading coefficient is det(G_off0)

theorem det_param_leading_coeff (N : ℕ) (L : ℝ) (K : ℕ) (z : Fin K → ℂ) :
    ∃ (coeffs : Fin (2 * N + 1) → ℝ),
      (∀ (w : ℝ), det_param N L K z w = Finset.sum Finset.univ
        (fun i : Fin (2 * N + 1) => coeffs i * w ^ (i : ℕ))) ∧
      coeffs ⟨2 * N, by omega⟩ = (G_off0 N L K z).det := by
  -- The leading term of det(wG - K) is w^(2N) * det(G) by multilinearity
  sorry

#print axioms det_param_poly_degree
#print axioms det_param_leading_coeff
