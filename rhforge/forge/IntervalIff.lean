import IntervalWeightedNonneg
import IntervalConverse

/-!
# §1811 — Interval localization as a single iff

For a real root multiset `R` with all roots simple, the `(4−r²)`-weighted Hankel form is positive semidefinite
(at every degree, every real vector) **iff** every root lies in `[−2,2]` (`r² ≤ 4`).  Combines §1810 (⟸),
§1809 (⟹), and the real decomposition `R = r₀ ::ₘ (R.erase r₀)`.  The interval analog of §1806; together they
give "all roots real AND in [−2,2] ⟺ both Hankel forms PSD", the exact input of the on-circle criterion §1789.
EXIT0 / std-axioms.
-/

open Finset

/-- **Interval iff.**  All-simple real `R`: `(4−r²)`-weighted Hankel form PSD ⟺ all roots in `[−2,2]`. -/
theorem interval_weighted_psd_iff (R : Multiset ℝ) (hsimple : ∀ r ∈ R, R.count r = 1) :
    (∀ (n : ℕ) (c : Fin (n + 1) → ℝ),
        0 ≤ (R.map (fun r => (4 - r ^ 2) * (∑ i : Fin (n + 1), c i * r ^ (i : ℕ)) ^ 2)).sum)
      ↔ (∀ r ∈ R, r ^ 2 ≤ 4) := by
  constructor
  · intro hpsd
    by_contra hnot
    push_neg at hnot
    obtain ⟨r₀, hr₀mem, hr₀out⟩ := hnot
    have hout : 4 < r₀ ^ 2 := hr₀out
    have hReq : R = r₀ ::ₘ R.erase r₀ := (Multiset.cons_erase hr₀mem).symm
    have hr₀rest : r₀ ∉ R.erase r₀ := by
      intro h
      have hp := Multiset.count_pos.mpr h
      rw [Multiset.count_erase_self, hsimple r₀ hr₀mem] at hp; omega
    obtain ⟨n, c, hval⟩ := interval_form_neg_of_outside (R.erase r₀) r₀ hr₀rest hout
    have hge := hpsd n c
    rw [hReq] at hge
    linarith
  · intro hin n c
    exact interval_weighted_nonneg R c hin
