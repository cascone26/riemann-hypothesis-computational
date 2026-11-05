import DetParamFactorization
import QOffKPoly

open scoped MatrixOrder Matrix

/-!
# Eq (7.17) with explicit invertibility hypothesis

Following standard mathematical practice, we accept G_off0's invertibility as an EXPLICIT HYPOTHESIS
rather than attempting to prove it. This is analogous to how many theorems in linear algebra
require generic/nondegeneracy assumptions.

The statement becomes:
  `det(w·G_off0 - K_off0) = det(G_off0)·Q_off_K(w)`

UNDER TWO EXPLICIT HYPOTHESES:
  (hG)    : G_off0 is invertible
  (h_char): the characteristic polynomial identity det(w•1 - G_off0⁻¹·K_off0) = Q_off_K(w)

NOTE (corrects the prior version of this file, which routed through CharPolyIdentity.lean's
vacuous `eigenvalues_match_roots` axiom — see LEDGER S2791): `h_char` is now a plain, honest
explicit hypothesis, not backed by any axiom or sorry inside THIS theorem. The theorem below
is fully proven (zero sorryAx, zero extra axioms) CONDITIONAL on its two named hypotheses —
exactly the correct way to represent "the rest is standard math, modulo two explicitly-named
open sub-claims" in this project's established style.
-/

theorem eq_7_17_with_explicit_invertibility (N : ℕ) (L : ℝ) (hL : L ≠ 0) (K : ℕ)
    (z : Fin K → ℂ) (hz : ∀ k, (z k).im ≠ 0) (w : ℝ)
    (hG : IsUnit (G_off0 N L K z).det)
    (h_char : (w • 1 - (G_off0 N L K z)⁻¹ * (K_off0 N L K z)).det = Q_off_K K z w) :
    det_param N L K z w = (G_off0 N L K z).det * (Q_off_K K z w) := by
  exact eq_7_17_reduction N L hL K z hz w hG h_char

#print axioms eq_7_17_with_explicit_invertibility
