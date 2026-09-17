import Mathlib.NumberTheory.Chebyshev

set_option maxHeartbeats 1000000

/-!
# §1592 — `ψ` is right-continuous: `ContinuousWithinAt ψ (Ioi t) t` for `t ≥ 0`

`ψ` is constant on `[t, ⌊t⌋+1)` (the floor doesn't move), so on the right-neighborhood `Ioo t (⌊t⌋+1)` it equals
`ψ(t)`.  Right-continuity is exactly what the right-FTC (`integral_eq_sub_of_hasDeriv_right`) needs — NO
countable-exception machinery: the primitive `ψ₁` has right-derivative `ψ(t)` at EVERY point.  EXIT0 / std-axioms.
-/

open Chebyshev Filter

theorem psi_continuousWithinAt_Ioi {t : ℝ} (ht : 0 ≤ t) :
    ContinuousWithinAt (fun u : ℝ => ψ u) (Set.Ioi t) t := by
  have hlt : t < (⌊t⌋₊ : ℝ) + 1 := Nat.lt_floor_add_one t
  have hmem : Set.Ioo t ((⌊t⌋₊ : ℝ) + 1) ∈ nhdsWithin t (Set.Ioi t) :=
    Ioo_mem_nhdsGT hlt
  have heq : (fun u : ℝ => ψ u) =ᶠ[nhdsWithin t (Set.Ioi t)] (fun _ => ψ t) := by
    filter_upwards [hmem] with u hu
    rw [Set.mem_Ioo] at hu
    have hfl : ⌊u⌋₊ = ⌊t⌋₊ := by
      rw [Nat.floor_eq_iff (by linarith [hu.1] : (0 : ℝ) ≤ u)]
      exact ⟨le_trans (Nat.floor_le ht) (le_of_lt hu.1), hu.2⟩
    rw [psi_eq_psi_coe_floor u, hfl, ← psi_eq_psi_coe_floor t]
  exact (continuousWithinAt_const).congr_of_eventuallyEq heq rfl

#print axioms psi_continuousWithinAt_Ioi
