import Mathlib.Analysis.SpecialFunctions.Integrals.Basic

set_option maxHeartbeats 1000000

/-!
# §1579 — strip integral identity: `(n+1)^{-s} - n^{-s} = -s·∫_n^{n+1} t^{-s-1} dt`

Second brick of the Mellin pillar.  From Mathlib's `integral_cpow` (`r = -s-1`, `r+1 = -s ≠ 0`,
`0 ∉ [[n, n+1]]` since `n ≥ 1`): `∫_n^{n+1} t^{-s-1} dt = ((n+1)^{-s} - n^{-s})/(-s)`.  Rearranged, the
`Δn^{-s}` differences of §1578 become integrals — the continuous half of the Mellin representation.
EXIT0 / std-axioms.
-/

open Complex

theorem cpow_diff_eq_integral {s : ℂ} (hs : s ≠ 0) {n : ℕ} (hn : 0 < n) :
    (((n + 1 : ℕ) : ℂ)) ^ (-s) - (n : ℂ) ^ (-s)
      = -s * ∫ t : ℝ in (n : ℝ)..(((n + 1 : ℕ) : ℝ)), (t : ℂ) ^ (-s - 1) := by
  have hnR : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
  have hr1 : (-s - 1 : ℂ) ≠ -1 := by
    intro h
    apply hs
    linear_combination -h
  have h0 : (0 : ℝ) ∉ Set.uIcc (n : ℝ) (((n + 1 : ℕ) : ℝ)) := by
    intro h
    rw [Set.mem_uIcc] at h
    push_cast at h
    rcases h with ⟨h1, _⟩ | ⟨h1, _⟩ <;> linarith
  rw [integral_cpow (Or.inr ⟨hr1, h0⟩)]
  have h1 : (-s - 1 : ℂ) + 1 = -s := by ring
  rw [h1]
  have hns : (-s : ℂ) ≠ 0 := neg_ne_zero.mpr hs
  push_cast
  field_simp

#print axioms cpow_diff_eq_integral
