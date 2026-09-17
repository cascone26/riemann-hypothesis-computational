import Mathlib.Algebra.BigOperators.Group.Multiset.Basic
import Mathlib.Data.Complex.Basic
import Mathlib.Tactic

/-!
# §1805 — Conjugation-closed decomposition around a simple non-real root

A conjugation-closed root multiset `R` with a SIMPLE (`count = 1`) non-real root `r₀` splits as
`R = r₀ ::ₘ conj r₀ ::ₘ rest` with `rest` conjugation-closed and `r₀, conj r₀ ∉ rest` — the data
§1803 consumes.  EXIT0 / std-axioms.
-/

open Multiset

private theorem conj_inj : Function.Injective (starRingEnd ℂ) :=
  Function.LeftInverse.injective (fun z => Complex.conj_conj z)

/-- **Conj-closed decomposition around a simple non-real root.** -/
theorem conj_closed_decomp (R : Multiset ℂ) (r₀ : ℂ)
    (hcc : R.map (starRingEnd ℂ) = R) (hmem : r₀ ∈ R) (him : r₀.im ≠ 0)
    (hsimple : R.count r₀ = 1) :
    ∃ rest : Multiset ℂ, R = r₀ ::ₘ (starRingEnd ℂ r₀) ::ₘ rest
      ∧ rest.map (starRingEnd ℂ) = rest
      ∧ r₀ ∉ rest ∧ (starRingEnd ℂ r₀) ∉ rest := by
  have hne : starRingEnd ℂ r₀ ≠ r₀ := by
    intro h; apply him
    have := congrArg Complex.im h
    simp only [Complex.conj_im] at this; linarith
  have hbarmem : starRingEnd ℂ r₀ ∈ R := by
    have : starRingEnd ℂ r₀ ∈ R.map (starRingEnd ℂ) := Multiset.mem_map_of_mem _ hmem
    rwa [hcc] at this
  have hbarmem' : starRingEnd ℂ r₀ ∈ R.erase r₀ :=
    Multiset.mem_erase_of_ne hne |>.mpr hbarmem
  have hcountbar : R.count (starRingEnd ℂ r₀) = 1 := by
    have h := Multiset.count_map_eq_count' (starRingEnd ℂ) R conj_inj r₀
    rw [hcc, hsimple] at h; exact h
  refine ⟨(R.erase r₀).erase (starRingEnd ℂ r₀), ?_, ?_, ?_, ?_⟩
  · rw [Multiset.cons_erase hbarmem', Multiset.cons_erase hmem]
  · rw [Multiset.map_erase _ conj_inj, Multiset.map_erase _ conj_inj, hcc, Complex.conj_conj,
        Multiset.erase_comm]
  · intro h
    have h1 := Multiset.count_pos.mpr (Multiset.mem_of_mem_erase h)
    rw [Multiset.count_erase_self, hsimple] at h1; omega
  · intro h
    have h1 := Multiset.count_pos.mpr h
    rw [Multiset.count_erase_self, Multiset.count_erase_of_ne hne, hcountbar] at h1; omega
