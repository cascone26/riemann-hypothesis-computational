import PNTPsiError
import Mathlib.NumberTheory.Chebyshev

set_option maxHeartbeats 1000000

open Chebyshev Real

/-!
# §1705 ★★★★ — THE PNT FOR θ: `|θ(x) − x| ≤ C·x·e^{-c√log x}`

`ψ − θ ≤ 2√x·log x` (Mathlib, Chebyshev.psi_sub_theta_le) is EXPONENTIALLY smaller than the dlVP
error scale: in `L = log x` form, `2L·e^{c√L} ≤ e^{L/4}·e^{L/4} = e^{L/2}` once `√L ≥ 4c` and
`L²/64 ≥ 2L` — so `2√x·log x ≤ x·e^{-c√log x}` for large `x` and the prime-power correction
disappears into the constant: `|θ(x) − x| ≤ (C+1)·x·e^{-c√log x}`.
-/

theorem two_L_exp_le {c' L : ℝ} (hc : 0 < c') (hL : max (16 * c' ^ 2) 128 ≤ L) :
    2 * Real.exp (L / 2) * L ≤ Real.exp L * Real.exp (-(c' * Real.sqrt L)) := by
  have hL128 : (128 : ℝ) ≤ L := le_trans (le_max_right _ _) hL
  have hL16 : 16 * c' ^ 2 ≤ L := le_trans (le_max_left _ _) hL
  have hL0 : (0 : ℝ) ≤ L := by linarith
  -- (i) c'·√L ≤ L/4
  have hsqL : 4 * c' ≤ Real.sqrt L := by
    have h1 : ((4 * c') : ℝ) ^ 2 ≤ L := by nlinarith
    have h2 := Real.sqrt_le_sqrt h1
    rwa [Real.sqrt_sq (by positivity)] at h2
  have hi : c' * Real.sqrt L ≤ L / 4 := by
    have h1 := mul_le_mul_of_nonneg_right hsqL (Real.sqrt_nonneg L)
    have h2 : Real.sqrt L * Real.sqrt L = L := Real.mul_self_sqrt hL0
    nlinarith [h1, h2]
  -- (ii) 2L ≤ e^{L/4}
  have hii : 2 * L ≤ Real.exp (L / 4) := by
    have h1 : 1 + L / 8 ≤ Real.exp (L / 8) := by
      have := Real.add_one_le_exp (L / 8)
      linarith
    have h2 : Real.exp (L / 4) = Real.exp (L / 8) ^ 2 := by
      rw [sq, ← Real.exp_add]; congr 1; ring
    nlinarith [h1, h2, Real.exp_pos (L / 8)]
  -- combine: 2L·e^{c'√L} ≤ e^{L/4}·e^{L/4} = e^{L/2}
  have hkey : 2 * L * Real.exp (c' * Real.sqrt L) ≤ Real.exp (L / 2) := by
    have h1 : Real.exp (c' * Real.sqrt L) ≤ Real.exp (L / 4) := Real.exp_le_exp.mpr hi
    have h2 : Real.exp (L / 2) = Real.exp (L / 4) * Real.exp (L / 4) := by
      rw [← Real.exp_add]; congr 1; ring
    rw [h2]
    exact mul_le_mul hii h1 (Real.exp_pos _).le (Real.exp_pos _).le
  -- multiply through by e^{L/2 - c'√L} > 0
  have hz : (0 : ℝ) < Real.exp (L / 2 + -(c' * Real.sqrt L)) := Real.exp_pos _
  have hmul := mul_le_mul_of_nonneg_right hkey hz.le
  calc 2 * Real.exp (L / 2) * L
      = 2 * L * Real.exp (c' * Real.sqrt L)
          * Real.exp (L / 2 + -(c' * Real.sqrt L)) := by
        rw [mul_assoc (2 * L), ← Real.exp_add]
        ring_nf
    _ ≤ Real.exp (L / 2) * Real.exp (L / 2 + -(c' * Real.sqrt L)) := hmul
    _ = Real.exp L * Real.exp (-(c' * Real.sqrt L)) := by
        rw [← Real.exp_add, ← Real.exp_add]
        congr 1
        ring

theorem pnt_theta_error :
    ∃ C c x₀ : ℝ, 0 < C ∧ 0 < c ∧ 1 ≤ x₀ ∧ ∀ x : ℝ, x₀ ≤ x →
      |θ x - x| ≤ C * x * Real.exp (-(c * Real.sqrt (Real.log x))) := by
  obtain ⟨C, c, x₀, hC, hc, hx₀, hbound⟩ := pnt_psi_error
  refine ⟨C + 1, c, max x₀ (Real.exp (max (16 * c ^ 2) 128)),
    by linarith, hc, le_trans hx₀ (le_max_left _ _), ?_⟩
  intro x hx
  have hxx₀ : x₀ ≤ x := le_trans (le_max_left _ _) hx
  have hxexp : Real.exp (max (16 * c ^ 2) 128) ≤ x := le_trans (le_max_right _ _) hx
  have hx0 : (0 : ℝ) < x := lt_of_lt_of_le (Real.exp_pos _) hxexp
  have hx1 : (1 : ℝ) ≤ x := by linarith [hx₀]
  have hψ := hbound x hxx₀
  have hL : max (16 * c ^ 2) 128 ≤ Real.log x := (Real.le_log_iff_exp_le hx0).mpr hxexp
  -- 2√x·log x ≤ x·e^{-c√log x}
  have hcross : 2 * Real.sqrt x * Real.log x
      ≤ x * Real.exp (-(c * Real.sqrt (Real.log x))) := by
    have h1 := two_L_exp_le hc hL
    have hsx : Real.exp (Real.log x / 2) = Real.sqrt x := by
      rw [Real.sqrt_eq_rpow, Real.rpow_def_of_pos hx0]
      congr 1
      ring
    have hex : Real.exp (Real.log x) = x := Real.exp_log hx0
    rw [hsx, hex] at h1
    exact h1
  -- θ between ψ − 2√x log x and ψ
  have hθψ : θ x ≤ ψ x := theta_le_psi x
  have hψθ : ψ x - θ x ≤ 2 * Real.sqrt x * Real.log x := psi_sub_theta_le hx1
  have hψpair := abs_le.mp hψ
  rw [abs_le]
  constructor
  · -- −(C+1)xE ≤ θ x − x : θ x ≥ ψ x − 2√x log x ≥ x − CxE − xE
    have hxE0 : (0 : ℝ) ≤ x * Real.exp (-(c * Real.sqrt (Real.log x))) := by positivity
    nlinarith [hψpair.1, hψθ, hcross]
  · -- θ x − x ≤ ψ x − x ≤ CxE ≤ (C+1)xE
    have hxE0 : (0 : ℝ) ≤ x * Real.exp (-(c * Real.sqrt (Real.log x))) := by positivity
    nlinarith [hψpair.2, hθψ]

#print axioms two_L_exp_le
#print axioms pnt_theta_error
