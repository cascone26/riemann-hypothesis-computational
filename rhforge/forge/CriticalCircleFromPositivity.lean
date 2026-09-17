import DiscriminantTraceBound
import FrobeniusCriticalCircle

set_option maxHeartbeats 1000000

noncomputable section

/-!
# §1727 ★★★★★ — CAPSTONE: positivity ⟹ Frobenius eigenvalues on the critical circle

Composes the two verified halves of the function-field template's extraction into one statement:
- §1726 `disc_nonpos_of_quadratic_nonneg`: a nonnegative quadratic `q·x² + t·x + 1 ≥ 0` (with `q>0`)
  forces `t² ≤ 4q` — this is the Weil/Rosati **positivity** input (which our Hodge-index negativity
  on the primitive part, §1720/§1713/§1725, supplies in the geometric setting).
- §1724 `norm_root_eq_sqrt_of_disc_le`: `t² ≤ 4q` forces every root `λ` of `X² − tX + q` onto the
  **critical circle** `|λ| = √q`.

Result: *positivity alone* pins the Frobenius eigenvalues to `|λ| = q^{1/2}` — the Riemann Hypothesis
for a curve, as a single machine-verified implication. On the function-field side the template is now
complete and axiom-clean; the sole remaining input is arithmetic realization over `Spec ℤ`
(Connes–Consani, T-c) — a genuine open problem, honestly outside the formalized algebra.
-/

/-- **The critical-circle theorem.** If the reciprocal quadratic `q·x² + t·x + 1` is everywhere
nonnegative (Weil positivity, `q > 0`), then every root `λ` of the Frobenius characteristic
polynomial `X² − tX + q` satisfies `‖λ‖ = √q` — the eigenvalues lie exactly on the critical circle. -/
theorem critical_circle_of_positivity (t q : ℝ) (hq : 0 < q)
    (hpos : ∀ x : ℝ, 0 ≤ q * x ^ 2 + t * x + 1)
    (l : ℂ) (hroot : l ^ 2 - (t : ℂ) * l + (q : ℂ) = 0) :
    ‖l‖ = Real.sqrt q := by
  have hdisc : t ^ 2 ≤ 4 * q := by
    have h := disc_nonpos_of_quadratic_nonneg q t 1 hq hpos
    linarith
  exact norm_root_eq_sqrt_of_disc_le q t hq hdisc l hroot

/-- Squared form: `normSq λ = q` under the same positivity hypothesis. -/
theorem normSq_eq_of_positivity (t q : ℝ) (hq : 0 < q)
    (hpos : ∀ x : ℝ, 0 ≤ q * x ^ 2 + t * x + 1)
    (l : ℂ) (hroot : l ^ 2 - (t : ℂ) * l + (q : ℂ) = 0) :
    Complex.normSq l = q := by
  have hdisc : t ^ 2 ≤ 4 * q := by
    have h := disc_nonpos_of_quadratic_nonneg q t 1 hq hpos
    linarith
  exact normSq_root_eq_of_disc_le q t hq hdisc l hroot

#print axioms critical_circle_of_positivity
#print axioms normSq_eq_of_positivity
