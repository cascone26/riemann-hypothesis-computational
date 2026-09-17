import IntervalConverseCore
import RealRootIsolation

/-!
# §1809 — Interval converse: a real root outside [−2,2] breaks the weighted form

Composes §1808 (isolate a real root) with §1807 (an isolated root outside `[−2,2]` drives the `(4−r²)`-weighted
Hankel form negative).  Companion to §1791 (`hankel_interval_form_nonneg`, the ⟸): for real `r₀ ∉ rest` with
`r₀² > 4`, there is a real vector making the weighted form over `r₀ ::ₘ rest` strictly negative — the form is
not PSD.  EXIT0 / std-axioms.
-/

open Finset

/-- **Interval converse.**  A real root outside `[−2,2]` makes the `(4−r²)`-weighted Hankel form non-PSD. -/
theorem interval_form_neg_of_outside (rest : Multiset ℝ) (r₀ : ℝ)
    (h : r₀ ∉ rest) (hout : 4 < r₀ ^ 2) :
    ∃ (n : ℕ) (c : Fin (n + 1) → ℝ),
      ((r₀ ::ₘ rest).map
        (fun r => (4 - r ^ 2) * (∑ i : Fin (n + 1), c i * r ^ (i : ℕ)) ^ 2)).sum < 0 := by
  obtain ⟨n, c, hr₀, hrest⟩ := real_root_isolation rest r₀ h
  exact ⟨n, c, interval_isolated_neg rest c r₀ hr₀ hrest hout⟩
