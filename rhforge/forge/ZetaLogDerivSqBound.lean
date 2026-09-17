import NormLogDerivZetaBall
import ZetaSphereBoundCenter
import ZetaLowerBound
import PSeriesUpperBound
import LandauZeroDistance

set_option maxHeartbeats 2000000

open Complex Metric

/-!
# §1687 ★★★★★ — `ζ'/ζ = O(log²T)` AT THE LANDAU BALL — the classical bound, fully machine-checked

THE assembled bound: there are `c₁, A > 0` such that for all `T ≥ 6` and every `z` within `1/64`
of the Landau center `c(T) = 1 + 1/log(T+2) + iT` with `Re z ≥ 1 - c₁/log(T+2)` (the zero-free
side), `‖ζ'/ζ(z)‖ ≤ A·log²(T+2)`. Every ingredient banked tonight: §1683 (the ball bound),
§1684 (`M₀ = 5T` sphere growth), §1678+§1685 (`‖ζ(c)‖ ≥ 1/(1+log(T+2))` via Möbius + p-series),
§1686 (`δ = c₁/log(T+2)` from the zero-free region). This is the de la Vallée Poussin
`ζ'/ζ = O(log²T)` — the last analytic input of PNT-with-classical-error — as an axiom-clean
Lean theorem about Mathlib's `riemannZeta`.
-/

theorem norm_logDeriv_zeta_sq_bound :
    ∃ c₁ A : ℝ, 0 < c₁ ∧ 0 < A ∧ ∀ T : ℝ, 6 ≤ T → ∀ z : ℂ,
      ‖z - (((1 + 1 / Real.log (T + 2) : ℝ) : ℂ) + (T : ℝ) * Complex.I)‖ ≤ 1 / 64 →
      1 - c₁ / Real.log (T + 2) ≤ z.re →
      ‖logDeriv riemannZeta z‖ ≤ A * (Real.log (T + 2)) ^ 2 := by
  obtain ⟨c₁, hc₁, hdist⟩ := landau_ball_zero_distance
  have hlog3 : (0 : ℝ) < Real.log 3 := Real.log_pos (by norm_num)
  refine ⟨c₁, 4 / (c₁ * Real.log 3) + 256, hc₁, by positivity, ?_⟩
  intro T hT z hz hzre
  set L : ℝ := Real.log (T + 2) with hLdef
  have hL2 : 2 ≤ L := by
    rw [hLdef]
    rw [Real.le_log_iff_exp_le (by linarith)]
    have he := Real.exp_one_lt_d9
    calc Real.exp 2 = Real.exp 1 * Real.exp 1 := by
          rw [← Real.exp_add]; norm_num
      _ ≤ 2.7182818286 * 2.7182818286 := by
          apply mul_le_mul he.le he.le (Real.exp_pos 1).le (by norm_num)
      _ ≤ T + 2 := by nlinarith
  have hLpos : 0 < L := by linarith
  set c : ℂ := ((1 + 1 / L : ℝ) : ℂ) + (T : ℝ) * Complex.I with hcdef
  have hcre : c.re = 1 + 1 / L := by rw [hcdef]; simp
  have hcim : c.im = T := by rw [hcdef]; simp
  have hcre1 : 1 < c.re := by
    rw [hcre]
    have : 0 < 1 / L := by positivity
    linarith
  -- ball avoids the pole
  have hball : (1 : ℂ) ∉ ball c (1 / 2) := by
    intro h
    rw [mem_ball] at h
    have him := abs_im_le_norm ((1 : ℂ) - c)
    rw [show ((1 : ℂ) - c).im = -c.im by simp, hcim, abs_neg, abs_of_pos (by linarith)] at him
    rw [dist_eq_norm] at h
    linarith
  -- sphere bound M₀ = 5T
  have hε : (0 : ℝ) < 1 / L := by positivity
  have hε' : 1 / L ≤ 1 / 2 := one_div_le_one_div_of_le (by norm_num) hL2
  have hζsphere : ∀ s ∈ sphere c (1 / 4), ‖riemannZeta s‖ ≤ 5 * T := by
    intro s hs
    exact zeta_sphere_bound_center (by linarith) hε hε' s hs
  -- z in the R/16-ball
  have hzball : z ∈ closedBall c (1 / 4 / 16) := by
    rw [mem_closedBall, dist_eq_norm]
    calc ‖z - c‖ ≤ 1 / 64 := hz
      _ = 1 / 4 / 16 := by norm_num
  -- the δ-input from §1686
  have hdist' : ∀ ρ ∈ closedBall c (1 / 4 / 4), riemannZeta ρ = 0 → c₁ / L ≤ ‖z - ρ‖ := by
    intro ρ hρ hρ0
    apply hdist T (by linarith) z ρ hzre _ hρ0
    have him := abs_im_le_norm (ρ - c)
    rw [mem_closedBall, dist_eq_norm] at hρ
    rw [show (ρ - c).im = ρ.im - c.im by simp, hcim] at him
    have h1 : |ρ.im - T| ≤ 1 / 4 / 4 := le_trans him hρ
    have := abs_le.mp h1
    rw [abs_le]
    constructor <;> linarith
  -- the ball bound §1683
  have hmain := norm_logDeriv_riemannZeta_le_of_ball (by norm_num : (0:ℝ) < 1/4)
    (by norm_num : (1:ℝ)/4 < 1/2) (by positivity : (0:ℝ) < c₁ / L) hcre1 hball
    hζsphere hzball hdist'
  -- lower bound for ‖ζ(c)‖
  have hζclower : (1 + L)⁻¹ ≤ ‖riemannZeta c‖ := by
    have h678 := inv_tsum_le_norm_riemannZeta (s := c) (by rw [hcre]; linarith [hε])
    have h685 : ∑' n : ℕ, 1 / (n : ℝ) ^ c.re ≤ 1 + L := by
      have h := tsum_one_div_nat_rpow_le (σ := c.re) (by rw [hcre]; linarith [hε])
      have he : (1 : ℝ) + 1 / (c.re - 1) = 1 + L := by
        rw [hcre, show (1 : ℝ) + 1 / L - 1 = 1 / L by ring, one_div_one_div]
      linarith
    have hsumpos : (0 : ℝ) < ∑' n : ℕ, 1 / (n : ℝ) ^ c.re := by
      have h1le : 1 / ((1 : ℕ) : ℝ) ^ c.re ≤ ∑' n : ℕ, 1 / (n : ℝ) ^ c.re :=
        (Real.summable_one_div_nat_rpow.mpr (by rw [hcre]; linarith [hε])).le_tsum 1
          (fun j _ => by positivity)
      simp only [Nat.cast_one, Real.one_rpow, div_one] at h1le
      linarith
    calc (1 + L)⁻¹ ≤ (∑' n : ℕ, 1 / (n : ℝ) ^ c.re)⁻¹ := by
          rw [inv_eq_one_div, inv_eq_one_div]
          exact one_div_le_one_div_of_le hsumpos h685
      _ ≤ ‖riemannZeta c‖ := h678
  have hζcpos : 0 < ‖riemannZeta c‖ :=
    lt_of_lt_of_le (by positivity) hζclower
  -- log((5T+1)/‖ζ c‖) ≤ 4L
  have hlogbound : Real.log ((5 * T + 1) / ‖riemannZeta c‖) ≤ 4 * L := by
    have h1 : (5 * T + 1) / ‖riemannZeta c‖ ≤ (5 * T + 1) * (1 + L) := by
      rw [div_le_iff₀ hζcpos]
      have hcancel : ((1 : ℝ) + L)⁻¹ * (1 + L) = 1 := inv_mul_cancel₀ (by positivity)
      calc (5 * T + 1) = (5 * T + 1) * (((1 : ℝ) + L)⁻¹ * (1 + L)) := by
            rw [hcancel, mul_one]
        _ = (5 * T + 1) * (1 + L) * ((1 : ℝ) + L)⁻¹ := by ring
        _ ≤ (5 * T + 1) * (1 + L) * ‖riemannZeta c‖ := by gcongr
    have h2 : Real.log ((5 * T + 1) / ‖riemannZeta c‖) ≤ Real.log ((5 * T + 1) * (1 + L)) :=
      Real.log_le_log (by positivity) h1
    have h3 : Real.log ((5 * T + 1) * (1 + L))
        = Real.log (5 * T + 1) + Real.log (1 + L) := by
      rw [Real.log_mul (by positivity) (by positivity)]
    have h4 : Real.log (5 * T + 1) ≤ 2 * L := by
      have ha : (5 * T + 1 : ℝ) ≤ (T + 2) ^ 2 := by nlinarith
      have hb : Real.log (5 * T + 1) ≤ Real.log ((T + 2) ^ 2) :=
        Real.log_le_log (by linarith) ha
      rw [Real.log_pow] at hb
      push_cast at hb
      rw [hLdef]
      linarith
    have h5 : Real.log (1 + L) ≤ L := by
      have ha : (1 + L : ℝ) ≤ Real.exp L := by
        have := Real.add_one_le_exp L
        linarith
      calc Real.log (1 + L) ≤ Real.log (Real.exp L) :=
            Real.log_le_log (by linarith) ha
        _ = L := Real.log_exp L
    linarith
  -- final arithmetic
  calc ‖logDeriv riemannZeta z‖
      ≤ Real.log ((5 * T + 1) / ‖riemannZeta c‖) * (1 / (c₁ / L * Real.log 3) + 32 / (1 / 4)) :=
        hmain
    _ ≤ 4 * L * (1 / (c₁ / L * Real.log 3) + 32 / (1 / 4)) := by
        apply mul_le_mul_of_nonneg_right hlogbound
        positivity
    _ = 4 * L * (L / (c₁ * Real.log 3) + 128) := by
        congr 2
        · rw [show (c₁ / L * Real.log 3) = c₁ * Real.log 3 / L from by ring, one_div_div]
        · norm_num
    _ = 4 / (c₁ * Real.log 3) * L ^ 2 + 512 * L := by
        ring
    _ ≤ (4 / (c₁ * Real.log 3) + 256) * L ^ 2 := by
        nlinarith [hL2, hLpos]

#print axioms norm_logDeriv_zeta_sq_bound
