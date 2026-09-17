import SumEulerMaclaurin
import ZetaStripRep
import ZetaFractIntegral

/-!
# §1286 — Truncated ζ representation `ζ(s) = ∑_{k≤N} k^{-s} + N^{1-s}/(s-1) - s∫_N^∞ {t}t^{-(s+1)}`

The culmination of the truncated-representation chain (§1280–§1285).  Reconciling the finite Euler–Maclaurin
identity §1285 with the forge whole-range representation `ζ(s) = s/(s-1) - s∫_1^∞{t}t^{-(s+1)}` — splitting
`∫_1^∞ = ∫_1^N + ∫_N^∞` — the closed-form main terms collapse (by hand: `M = N^{1-s}/(s-1)`) to the standard
truncated Euler–Maclaurin representation, valid for `Re s > 0`, `s ≠ 1`:

  `ζ(s) = (∑_{k=2}^N k^{-s}) + 1 + N^{1-s}/(s-1) - s∫_N^∞ {t}t^{-(s+1)}`.

With `N ≍ |t|` the main-sum + `O(1)` remainder gives the log-growth bound §1269 → sharp zero-free region.
EXIT0 / std-axioms.
-/

open Finset MeasureTheory Set

theorem riemannZeta_eq_partial_sum_add_tail {N : ℕ} (hN : 1 ≤ N) {s : ℂ}
    (hs0 : 0 < s.re) (hs1 : s ≠ 1) :
    riemannZeta s = (∑ k ∈ Ioc 1 N, (k : ℂ) ^ (-s)) + 1 + (N : ℂ) ^ (1 - s) / (s - 1)
      - s * ∫ t in Ioi (N : ℝ), ((Int.fract t : ℝ) : ℂ) * (t : ℂ) ^ (-(s + 1)) := by
  have hs0' : s ≠ 0 := by
    intro h; rw [h] at hs0; simp at hs0
  have hsm1 : s - 1 ≠ 0 := sub_ne_zero.mpr hs1
  have h1ms : (1 : ℂ) - s ≠ 0 := fun h => hsm1 (by linear_combination -h)
  have hNR : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast hN
  have hrep := riemannZeta_eq_sub_fract_integral_strip hs0 hs1
  have hint := integrableOn_fract_cpow hs0
  have hsplit : (∫ t in Ioi (1 : ℝ), ((Int.fract t : ℝ) : ℂ) * (t : ℂ) ^ (-(s + 1)))
      = (∫ t in Ioc (1 : ℝ) (N : ℝ), ((Int.fract t : ℝ) : ℂ) * (t : ℂ) ^ (-(s + 1)))
        + (∫ t in Ioi (N : ℝ), ((Int.fract t : ℝ) : ℂ) * (t : ℂ) ^ (-(s + 1))) := by
    rw [← Set.Ioc_union_Ioi_eq_Ioi hNR]
    exact setIntegral_union (Set.Ioc_disjoint_Ioi le_rfl) measurableSet_Ioi
      (hint.mono_set Set.Ioc_subset_Ioi_self) (hint.mono_set (Set.Ioi_subset_Ioi hNR))
  have hem := sum_cpow_neg_euler_maclaurin hN hs0' hs1
  simp only [show (-s - 1 : ℂ) = -(s + 1) from by ring] at hem
  have hNc : (N : ℂ) ≠ 0 := by exact_mod_cast Nat.one_le_iff_ne_zero.mp hN
  have hNpow : (N : ℂ) ^ (1 - s) = (N : ℂ) ^ (-s) * (N : ℂ) := by
    rw [show (1 - s : ℂ) = -s + 1 from by ring, Complex.cpow_add _ _ hNc, Complex.cpow_one]
  rw [hrep, hsplit, hem, hNpow]
  field_simp
  ring

#print axioms riemannZeta_eq_partial_sum_add_tail
