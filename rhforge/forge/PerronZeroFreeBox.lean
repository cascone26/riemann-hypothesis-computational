import ZetaZeroFreeUniform
import Mathlib.Analysis.SpecialFunctions.Log.Basic

set_option maxHeartbeats 1000000

open Complex Set

open scoped Interval

/-!
# §1667 — THE ZERO-FREE BOX: DLVP region ⟹ §1665/§1666's rectangle hypothesis, uniformly in `T`

Bridges the zero-free tower to the winding arc: there is a single constant `c > 0` such that for
EVERY `T ≥ 2` and every `σ₀ > 1`, the rectangle `[1 - c/log(T+2), σ₀] × [-T, T]` satisfies the
full §1665 hypothesis — no `ζ`-zeros (by the uniform DLVP region §1577: any zero at height `γ`,
`|γ| ≤ T`, has `β ≤ 1 - c₀/log(|γ|+2) < 1 - c/log(T+2) = ` the box's left edge, since
`2c ≤ c₀`), and the box misses `0` and `-1` (left edge stays `≥ ½` by `c ≤ log 2`). The Perron
shift §1666 can now run on the `T`-dependent box `σ₁(T) = 1 - c/log(T+2)` — the exact geometry
of PNT-with-classical-error.
-/

theorem perron_zero_free_box :
    ∃ c : ℝ, 0 < c ∧ ∀ T σ₀ : ℝ, 2 ≤ T → 1 < σ₀ →
      ∀ u ∈ [[1 - c / Real.log (T + 2), σ₀]] ×ℂ [[-T, T]],
        (u ≠ 1 → riemannZeta u ≠ 0) ∧ u ≠ 0 ∧ u ≠ -1 := by
  obtain ⟨c₀, hc₀, huni⟩ := riemannZeta_zero_free_region_uniform
  refine ⟨min (c₀ / 2) (Real.log 2), lt_min (by linarith) (Real.log_pos (by norm_num)), ?_⟩
  set c : ℝ := min (c₀ / 2) (Real.log 2) with hcdef
  intro T σ₀ hT hσ₀ u hu
  have hc : 0 < c := lt_min (by linarith) (Real.log_pos (by norm_num))
  have hchalf : c ≤ c₀ / 2 := min_le_left _ _
  have hclog : c ≤ Real.log 2 := min_le_right _ _
  have hL : 0 < Real.log (T + 2) := Real.log_pos (by linarith)
  have hlog4 : 2 * Real.log 2 ≤ Real.log (T + 2) := by
    have h4 : Real.log 4 ≤ Real.log (T + 2) := Real.log_le_log (by norm_num) (by linarith)
    have h42 : Real.log 4 = 2 * Real.log 2 := by
      rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.log_pow]
      push_cast
      ring
    linarith
  -- the left edge stays ≥ ½
  have hhalf : (1 : ℝ) / 2 ≤ 1 - c / Real.log (T + 2) := by
    have hkey : c / Real.log (T + 2) ≤ 1 / 2 := by
      rw [div_le_iff₀ hL]
      linarith
    linarith
  have hσ₁lt : 1 - c / Real.log (T + 2) < 1 := by
    have : 0 < c / Real.log (T + 2) := div_pos hc hL
    linarith
  -- unpack membership
  rw [Complex.mem_reProdIm,
    Set.uIcc_of_le (show (1 : ℝ) - c / Real.log (T + 2) ≤ σ₀ by linarith),
    Set.uIcc_of_le (show -T ≤ T by linarith)] at hu
  obtain ⟨hre, him⟩ := hu
  rw [Set.mem_Icc] at hre him
  refine ⟨fun _ => ?_, fun h => ?_, fun h => ?_⟩
  · -- no zeros in the box, via the uniform region
    intro hz
    have hz' : riemannZeta ((u.re : ℝ) + Complex.I * (u.im : ℝ)) = 0 := by
      rw [show ((u.re : ℝ) : ℂ) + Complex.I * (u.im : ℝ) = u by
        rw [mul_comm]; exact re_add_im u]
      exact hz
    have hβ := huni u.re u.im hz'
    have habs : |u.im| ≤ T := abs_le.mpr ⟨him.1, him.2⟩
    have hL2 : 0 < Real.log (|u.im| + 2) := Real.log_pos (by linarith [abs_nonneg u.im])
    have hd1 : c₀ / Real.log (T + 2) ≤ c₀ / Real.log (|u.im| + 2) :=
      div_le_div_of_nonneg_left hc₀.le hL2 (Real.log_le_log (by positivity) (by linarith))
    have hd3 : c / Real.log (T + 2) + c / Real.log (T + 2) ≤ c₀ / Real.log (T + 2) := by
      rw [← add_div]
      gcongr
      linarith
    have hpos : 0 < c / Real.log (T + 2) := div_pos hc hL
    linarith [hre.1]
  · have h' := congrArg Complex.re h
    simp at h'
    linarith [hre.1]
  · have h' := congrArg Complex.re h
    simp at h'
    linarith [hre.1]

#print axioms perron_zero_free_box
