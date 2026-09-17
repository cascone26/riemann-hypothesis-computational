import HankelFullConverse
import ComplexHankelNonnegReal
import ConjClosedDecomp

/-!
# §1806 — The Hermite criterion as a single iff (CAPSTONE)

For a conjugation-closed root multiset `R` with all roots simple, the complex Hankel form is positive
semidefinite (nonneg real part for every real coefficient vector, at every degree) **iff** all roots are real.
Combines §1804 (all-real ⟹ PSD), §1805 (decomposition around a simple non-real root), and §1803 (a non-real
root ⟹ the form hits `-2`).  This is the classical Hermite real-rootedness criterion, unconditional and
axiom-clean.  EXIT0 / std-axioms.
-/

open Finset

/-- **Hermite criterion (iff).**  Conj-closed, all-simple `R`: complex Hankel form PSD ⟺ all roots real. -/
theorem hankel_psd_iff_all_real (R : Multiset ℂ) (hcc : R.map (starRingEnd ℂ) = R)
    (hsimple : ∀ r ∈ R, R.count r = 1) :
    (∀ (n : ℕ) (c : Fin (n + 1) → ℝ),
        0 ≤ ((R.map (fun r => (∑ i : Fin (n + 1), (c i : ℂ) * r ^ (i : ℕ)) ^ 2)).sum).re)
      ↔ (∀ r ∈ R, r.im = 0) := by
  constructor
  · intro hpsd
    by_contra hnot
    push_neg at hnot
    obtain ⟨r₀, hr₀mem, hr₀im⟩ := hnot
    obtain ⟨rest, hReq, hrestcc, hr₀rest, hbarrest⟩ :=
      conj_closed_decomp R r₀ hcc hr₀mem hr₀im (hsimple r₀ hr₀mem)
    obtain ⟨n, c, hval⟩ := hankel_not_psd_of_nonreal rest r₀ hr₀im hrestcc hr₀rest hbarrest
    have hge := hpsd n c
    rw [hReq, hval, show (-2 : ℂ) = ((-2 : ℝ) : ℂ) by norm_cast, Complex.ofReal_re] at hge
    linarith
  · intro hreal n c
    exact complex_hankel_nonneg_of_real R c hreal
