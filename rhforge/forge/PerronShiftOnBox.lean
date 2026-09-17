import PerronVlineShift
import PerronZeroFreeBox

set_option maxHeartbeats 4000000

open Complex MeasureTheory Set

open scoped Interval

/-!
# §1668 ★ — THE PERRON SHIFT ON THE ZERO-FREE BOX: §1666 ∘ §1667, one theorem

The composition anchor for the PNT-error assembly: a single constant `c > 0` such that for every
`T ≥ 2`, `σ₀ > 1`, `x > 0`, the truncated Perron integral of `F` on `Re = σ₀` differs from the one
on the zero-free abscissa `Re = σ₁(T) = 1 - c/log(T+2)` by exactly `2π·(x²/2)` plus the two
horizontal connectors at heights `±T`. Everything downstream (the effective error term) is now
"bound the connectors on `[σ₁(T), σ₀] × {±T}`" — for which Mathlib's `borelCaratheodory` +
a Landau-lemma build is the loaded route.
-/

theorem perron_F_shift_on_box :
    ∃ c : ℝ, 0 < c ∧ ∀ T σ₀ x : ℝ, 2 ≤ T → 1 < σ₀ → 0 < x →
      (∫ t in (-T)..T,
          -logDeriv riemannZeta ((σ₀ : ℂ) + (t : ℝ) * I) *
              (x : ℂ) ^ (((σ₀ : ℂ) + (t : ℝ) * I) + 1) /
            (((σ₀ : ℂ) + (t : ℝ) * I) * (((σ₀ : ℂ) + (t : ℝ) * I) + 1))) -
        (∫ t in (-T)..T,
          -logDeriv riemannZeta (((1 - c / Real.log (T + 2) : ℝ) : ℂ) + (t : ℝ) * I) *
              (x : ℂ) ^ ((((1 - c / Real.log (T + 2) : ℝ) : ℂ) + (t : ℝ) * I) + 1) /
            ((((1 - c / Real.log (T + 2) : ℝ) : ℂ) + (t : ℝ) * I) *
              ((((1 - c / Real.log (T + 2) : ℝ) : ℂ) + (t : ℝ) * I) + 1)))
        = 2 * (Real.pi : ℂ) * ((x : ℂ) ^ (2 : ℂ) / 2) +
          I * ((∫ t in (1 - c / Real.log (T + 2))..σ₀,
              -logDeriv riemannZeta ((t : ℂ) + ((-T : ℝ) : ℂ) * I) *
                  (x : ℂ) ^ (((t : ℂ) + ((-T : ℝ) : ℂ) * I) + 1) /
                (((t : ℂ) + ((-T : ℝ) : ℂ) * I) * (((t : ℂ) + ((-T : ℝ) : ℂ) * I) + 1))) -
            (∫ t in (1 - c / Real.log (T + 2))..σ₀,
              -logDeriv riemannZeta ((t : ℂ) + ((T : ℝ) : ℂ) * I) *
                  (x : ℂ) ^ (((t : ℂ) + ((T : ℝ) : ℂ) * I) + 1) /
                (((t : ℂ) + ((T : ℝ) : ℂ) * I) * (((t : ℂ) + ((T : ℝ) : ℂ) * I) + 1)))) := by
  obtain ⟨c, hc, hbox⟩ := perron_zero_free_box
  refine ⟨c, hc, ?_⟩
  intro T σ₀ x hT hσ₀ hx
  set σ₁ : ℝ := 1 - c / Real.log (T + 2) with hσ₁def
  have hL : 0 < Real.log (T + 2) := Real.log_pos (by linarith)
  have hσ₁lt : σ₁ < 1 := by
    have : 0 < c / Real.log (T + 2) := div_pos hc hL
    rw [hσ₁def]
    linarith
  set z : ℂ := (σ₁ : ℂ) + ((-T : ℝ) : ℂ) * I with hzdef
  set w : ℂ := (σ₀ : ℂ) + ((T : ℝ) : ℂ) * I with hwdef
  have hz_re : z.re = σ₁ := by rw [hzdef]; simp
  have hz_im : z.im = -T := by rw [hzdef]; simp
  have hw_re : w.re = σ₀ := by rw [hwdef]; simp
  have hw_im : w.im = T := by rw [hwdef]; simp
  have h := perron_F_vline_shift hx z w
    (by rw [hz_re]; exact hσ₁lt) (by rw [hw_re]; exact hσ₀)
    (by rw [hz_im]; linarith) (by rw [hw_im]; linarith)
    (fun u hu => hbox T σ₀ hT hσ₀ u
      (by rw [hz_re, hw_re, hz_im, hw_im] at hu; exact hu))
  rw [hz_re, hz_im, hw_re, hw_im] at h
  exact h

#print axioms perron_F_shift_on_box
