import ZetaZeroFreeUniform

set_option maxHeartbeats 1000000

open Complex

/-!
# §1686 — the δ-DISCHARGE: zero-distance `c₁/log(T+2)` at the Landau ball

Discharges §1683's `hdist` hypothesis from the uniform DLVP zero-free region (§1577): there is
`c₁ > 0` such that for `T ≥ 2`, any point `z` with `Re z ≥ 1 - c₁/log(T+2)` keeps distance
`≥ c₁/log(T+2)` from EVERY ζ-zero of height `|Im ρ| ≤ T + 1`. Mechanism: such a zero has
`β ≤ 1 - c₀/log(|γ|+2) ≤ 1 - c₀/(2·log(T+2))` (since `|γ|+2 ≤ 2T+2 ≤ (T+2)²`-free doubling), so
with `c₁ = c₀/4` the horizontal gap alone is `≥ c₁/log(T+2)`. Purely real-arithmetic on top of the
zero-free tower.
-/

theorem landau_ball_zero_distance :
    ∃ c₁ : ℝ, 0 < c₁ ∧ ∀ T : ℝ, 2 ≤ T → ∀ z ρ : ℂ,
      1 - c₁ / Real.log (T + 2) ≤ z.re → |ρ.im| ≤ T + 1 → riemannZeta ρ = 0 →
      c₁ / Real.log (T + 2) ≤ ‖z - ρ‖ := by
  obtain ⟨c₀, hc₀, huni⟩ := riemannZeta_zero_free_region_uniform
  refine ⟨c₀ / 4, by linarith, ?_⟩
  intro T hT z ρ hzre hρim hρzero
  have hL : 0 < Real.log (T + 2) := Real.log_pos (by linarith)
  -- the zero-free bound at ρ
  have hz' : riemannZeta ((ρ.re : ℝ) + Complex.I * (ρ.im : ℝ)) = 0 := by
    rw [show ((ρ.re : ℝ) : ℂ) + Complex.I * (ρ.im : ℝ) = ρ by
      rw [mul_comm]; exact re_add_im ρ]
    exact hρzero
  have hβ := huni ρ.re ρ.im hz'
  -- log(|γ|+2) ≤ 2·log(T+2)
  have hlog2 : Real.log (|ρ.im| + 2) ≤ 2 * Real.log (T + 2) := by
    have h1 : |ρ.im| + 2 ≤ 2 * T + 2 := by linarith
    have h2 : Real.log (|ρ.im| + 2) ≤ Real.log (2 * T + 2) :=
      Real.log_le_log (by positivity) h1
    have h3 : Real.log (2 * T + 2) = Real.log 2 + Real.log (T + 1) := by
      rw [show 2 * T + 2 = 2 * (T + 1) by ring,
        Real.log_mul (by norm_num) (by linarith)]
    have h4 : Real.log 2 ≤ Real.log (T + 2) := Real.log_le_log (by norm_num) (by linarith)
    have h5 : Real.log (T + 1) ≤ Real.log (T + 2) := Real.log_le_log (by linarith) (by linarith)
    linarith
  -- hence β stays left of the line by c₀/(2·log(T+2))
  have hL2 : 0 < Real.log (|ρ.im| + 2) := Real.log_pos (by
    have := abs_nonneg ρ.im
    linarith)
  have hβ2 : ρ.re ≤ 1 - c₀ / (2 * Real.log (T + 2)) := by
    have hd : c₀ / (2 * Real.log (T + 2)) ≤ c₀ / Real.log (|ρ.im| + 2) :=
      div_le_div_of_nonneg_left hc₀.le hL2 hlog2
    linarith
  -- the horizontal gap
  have hgap : c₀ / 4 / Real.log (T + 2) ≤ z.re - ρ.re := by
    have e1 : c₀ / (2 * Real.log (T + 2)) = c₀ / 2 / Real.log (T + 2) := by
      ring
    rw [e1] at hβ2
    have e2 : c₀ / 2 / Real.log (T + 2) - c₀ / 4 / Real.log (T + 2)
        = c₀ / 4 / Real.log (T + 2) := by
      ring
    linarith
  calc c₀ / 4 / Real.log (T + 2) ≤ z.re - ρ.re := hgap
    _ ≤ |z.re - ρ.re| := le_abs_self _
    _ = |(z - ρ).re| := by simp
    _ ≤ ‖z - ρ‖ := abs_re_le_norm _

#print axioms landau_ball_zero_distance
