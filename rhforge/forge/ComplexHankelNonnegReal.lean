import ComplexHankel

/-!
# §1804 — Complex Hankel form is PSD when all roots are real (the ⟸ of the criterion)

If every root in `R` is real (`r.im = 0`), then for real coefficients `c` the complex Hankel form (§1795)
`∑_{r∈R}(∑_i c_i r^i)²` has nonnegative real part: a real root + real coeffs make each `∑_i c_i r^i` a real `W`,
so its square contributes `W² ≥ 0`.  With §1803 (`HankelFullConverse`) this is the full Hermite equivalence
PSD ⟺ all-roots-real.  EXIT0 / std-axioms.
-/

open Finset

variable {d : ℕ}

/-- **All-real ⟹ complex Hankel form PSD.** -/
theorem complex_hankel_nonneg_of_real (R : Multiset ℂ) (c : Fin (d + 1) → ℝ)
    (hR : ∀ r ∈ R, r.im = 0) :
    0 ≤ ((R.map (fun r => (∑ i : Fin (d + 1), (c i : ℂ) * r ^ (i : ℕ)) ^ 2)).sum).re := by
  induction R using Multiset.induction with
  | empty => simp
  | cons a R ih =>
    simp only [Multiset.map_cons, Multiset.sum_cons, Complex.add_re]
    have ha : a.im = 0 := hR a (Multiset.mem_cons_self a R)
    have hrest : ∀ r ∈ R, r.im = 0 := fun r hr => hR r (Multiset.mem_cons_of_mem hr)
    have ha' : a = ((a.re : ℝ) : ℂ) := by
      have h := Complex.re_add_im a
      rw [ha] at h; simp at h; exact h.symm
    have hapow : ∀ i : ℕ, a ^ i = ((a.re ^ i : ℝ) : ℂ) := by
      intro i; rw [Complex.ofReal_pow, ← ha']
    have hwW : (∑ i : Fin (d + 1), (c i : ℂ) * a ^ (i : ℕ))
        = (((∑ i : Fin (d + 1), c i * a.re ^ (i : ℕ)) : ℝ) : ℂ) := by
      rw [Complex.ofReal_sum]
      refine Finset.sum_congr rfl (fun i _ => ?_)
      rw [hapow, Complex.ofReal_mul]
    have hterm : 0 ≤ ((∑ i : Fin (d + 1), (c i : ℂ) * a ^ (i : ℕ)) ^ 2).re := by
      rw [hwW, ← Complex.ofReal_pow, Complex.ofReal_re]
      exact sq_nonneg _
    exact add_nonneg hterm (ih hrest)
