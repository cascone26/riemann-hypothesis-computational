import Mathlib.Analysis.Complex.Norm
import Mathlib.Analysis.SpecialFunctions.Pow.Real

/-!
# §1297 — Sharp bound on §1291's tail term keeping `‖s‖·N^{-σ}` paired

The §1296 caveat: bounding `N^{-σ} ≤ 1` in isolation loses the decay that cancels `‖s‖ ≈ |t|`, leaving a spurious
`|t|`-growth.  Keeping the product paired and using `|Im s| ≤ N` (which holds at `N = ⌈|t|⌉`), `‖s‖ ≤ σ + |Im s|`:

  `‖s‖·N^{-σ} ≤ σ + N^{1-σ}`     hence     `‖s‖·N^{-σ}/σ ≤ 1 + N^{1-σ}/σ`.

Now `N^{1-σ} ≤ e` in the regime `(1-σ)log N ≤ 1` (§1292) makes this genuinely bounded — no `|t|`-growth.  This is
the corrected tail bound feeding the sharp `‖ζ(σ+it)‖ ≲ log|t|`.  EXIT0 / std-axioms.
-/

theorem norm_mul_rpow_neg_div_le {s : ℂ} {N : ℝ} (hs0 : 0 < s.re) (hN : 1 ≤ N)
    (hsim : |s.im| ≤ N) :
    ‖s‖ * N ^ (-s.re) / s.re ≤ 1 + N ^ (1 - s.re) / s.re := by
  have hNpos : (0 : ℝ) < N := lt_of_lt_of_le one_pos hN
  have hpow_pos : (0 : ℝ) < N ^ (-s.re) := Real.rpow_pos_of_pos hNpos _
  have hpow_le1 : N ^ (-s.re) ≤ 1 := Real.rpow_le_one_of_one_le_of_nonpos hN (neg_nonpos.mpr hs0.le)
  have hsn : ‖s‖ ≤ s.re + |s.im| := by
    calc ‖s‖ ≤ |s.re| + |s.im| := Complex.norm_le_abs_re_add_abs_im s
      _ = s.re + |s.im| := by rw [abs_of_pos hs0]
  have hmul : N * N ^ (-s.re) = N ^ (1 - s.re) := by
    rw [show (1 - s.re : ℝ) = 1 + (-s.re) from by ring, Real.rpow_add hNpos, Real.rpow_one]
  have hkey : ‖s‖ * N ^ (-s.re) ≤ s.re + N ^ (1 - s.re) := by
    calc ‖s‖ * N ^ (-s.re) ≤ (s.re + |s.im|) * N ^ (-s.re) :=
          mul_le_mul_of_nonneg_right hsn hpow_pos.le
      _ = s.re * N ^ (-s.re) + |s.im| * N ^ (-s.re) := by ring
      _ ≤ s.re * 1 + N * N ^ (-s.re) := by
          apply add_le_add
          · exact mul_le_mul_of_nonneg_left hpow_le1 hs0.le
          · exact mul_le_mul_of_nonneg_right hsim hpow_pos.le
      _ = s.re + N ^ (1 - s.re) := by rw [mul_one, hmul]
  have hdiv : ‖s‖ * N ^ (-s.re) / s.re ≤ (s.re + N ^ (1 - s.re)) / s.re :=
    div_le_div_of_nonneg_right hkey hs0.le
  rwa [add_div, div_self (ne_of_gt hs0)] at hdiv

#print axioms norm_mul_rpow_neg_div_le
