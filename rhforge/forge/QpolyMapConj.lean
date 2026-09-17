import QpolyDef
import DicksonMapConj

/-!
# §1816 — `Qpoly` has real coefficients when the palindromic data is real

If every `v i` is real (`(v i).im = 0`), then `Polynomial.map conj (Qpoly v m) = Qpoly v m`: the coefficients
are real combinations of Dickson polynomials, which have real (conj-fixed, §1815) coefficients.  Hence
`(Qpoly v m).roots` is conjugation-closed — exactly the hypothesis §1806 (`hankel_psd_iff_all_real`) requires,
letting the on-circle criterion connect to the Hankel PSD iff.  EXIT0 / std-axioms.
-/

open Polynomial

/-- **`Qpoly` is real-coefficiented for real `v`.** -/
theorem Qpoly_map_conj (v : ℕ → ℂ) (m : ℕ) (hv : ∀ i, (v i).im = 0) :
    Polynomial.map (starRingEnd ℂ) (Qpoly v m) = Qpoly v m := by
  have hvc : ∀ i, starRingEnd ℂ (v i) = v i := by
    intro i; apply Complex.ext
    · simp
    · simp [hv i]
  simp only [Qpoly, Polynomial.map_add, Polynomial.map_C, Polynomial.map_sum,
    Polynomial.map_mul, dickson_map_conj, hvc]
