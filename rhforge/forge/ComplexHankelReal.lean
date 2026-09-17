import ComplexHankel

/-!
# §1796 — Complex Hankel form is REAL for a conjugation-closed root set

Converse assembly, step 1.  If the complex root multiset `R` is closed under conjugation (`R.map conj = R`,
true for the full root set of a REAL polynomial) then, for real coefficients `c`, the value of the complex
Hankel form `∑_{r∈R} (P_c r)²` is real (`im = 0`) — because complex conjugation permutes `R` and fixes real
`c`, so the form equals its own conjugate.  This lets "PSD of the complex Hankel form" be a meaningful (real)
hypothesis, which §1794 (conj-pair negativity) then refutes whenever a non-real root is present.
EXIT0 / std-axioms.
-/

open Finset

variable {d : ℕ}

/-- **Complex Hankel form is real** for a conjugation-closed root multiset and real coefficients. -/
theorem complex_hankel_form_real (R : Multiset ℂ) (c : Fin (d + 1) → ℝ)
    (hR : R.map (starRingEnd ℂ) = R) :
    ((R.map (fun r => (∑ i : Fin (d + 1), (c i : ℂ) * r ^ (i : ℕ)) ^ 2)).sum).im = 0 := by
  set P := fun r : ℂ => ∑ i : Fin (d + 1), (c i : ℂ) * r ^ (i : ℕ) with hP
  have hconj : ∀ r, starRingEnd ℂ (P r) = P (starRingEnd ℂ r) := by
    intro r
    simp only [hP, map_sum, map_mul, map_pow, Complex.conj_ofReal]
  have hSS : starRingEnd ℂ ((R.map (fun r => (P r) ^ 2)).sum)
      = (R.map (fun r => (P r) ^ 2)).sum := by
    rw [map_multiset_sum, Multiset.map_map]
    conv_rhs => rw [← hR, Multiset.map_map]
    refine congrArg Multiset.sum (Multiset.map_congr rfl (fun r _ => ?_))
    simp only [Function.comp_apply, map_pow, hconj]
  have himg := congrArg Complex.im hSS
  simp only [Complex.conj_im] at himg
  linarith
