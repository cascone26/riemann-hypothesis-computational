import Mathlib
import GammaVerticalDecay

/-!
# S2975 — exact modulus of ζ on the imaginary axis (lane 1: convexity endpoint)

From Mathlib's functional equation `ζ(1−s) = 2(2π)^{−s} Γ(s) cos(πs/2) ζ(s)` at `s = 1 − it`, the reflection-based
forge identity `‖Γ(it)‖² = π/(t sinh πt)` (GammaVerticalDecay), and `cos(π/2 − ix) = i sinh x`:
* `norm_Gamma_one_sub_mul_I_sq` — `‖Γ(1 − it)‖² = πt / sinh(πt)`  (t > 0).
* `norm_zeta_mul_I_sq` — `‖ζ(it)‖² · (2π cosh(πt/2)) = t · sinh(πt/2) · ‖ζ(1 − it)‖²`,
  i.e. `|ζ(it)|² = (t/2π)·tanh(πt/2)·|ζ(1−it)|²` EXACTLY.
* `norm_zeta_mul_I_sq_le` — `‖ζ(it)‖² ≤ (t/2π)·‖ζ(1 − it)‖²`.
This is the Re s = 0 edge of the convexity strip: with |ζ(1+it)| ≪ log t it gives |ζ(it)| ≪ |t|^{1/2} log t
(Phragmén–Lindelöf then interpolates to |ζ(σ+it)| ≪ |t|^{(1−σ)/2} log t — the next lane-1 brick).
-/

open Complex

namespace ZetaImagAxis

theorem norm_Gamma_one_sub_mul_I_sq (t : ℝ) (ht : 0 < t) :
    ‖Complex.Gamma (1 - (t : ℂ) * I)‖ ^ 2 = Real.pi * t / Real.sinh (Real.pi * t) := by
  have hs0 : ((t : ℂ) * I) ≠ 0 := mul_ne_zero (by exact_mod_cast ht.ne') I_ne_zero
  have h1s : Complex.Gamma (1 - (t : ℂ) * I)
      = (-((t : ℂ) * I)) * Complex.Gamma (-((t : ℂ) * I)) := by
    rw [show (1 : ℂ) - (t : ℂ) * I = (-((t : ℂ) * I)) + 1 by ring,
      Complex.Gamma_add_one (-((t : ℂ) * I)) (neg_ne_zero.mpr hs0)]
  have hconj : (starRingEnd ℂ) ((t : ℂ) * I) = -((t : ℂ) * I) := by
    simp [map_mul, Complex.conj_ofReal, Complex.conj_I]
  have hgconj : Complex.Gamma (-((t : ℂ) * I)) = (starRingEnd ℂ) (Complex.Gamma ((t : ℂ) * I)) := by
    rw [← hconj, Complex.Gamma_conj]
  have hsinhpos : 0 < Real.sinh (Real.pi * t) := Real.sinh_pos_iff.mpr (by positivity)
  rw [h1s, hgconj, norm_mul, norm_neg, Complex.norm_conj, norm_mul, Complex.norm_real,
    Complex.norm_I, mul_one, Real.norm_of_nonneg ht.le, mul_pow, norm_Gamma_mul_I_sq t ht]
  field_simp

theorem norm_zeta_mul_I_sq (t : ℝ) (ht : 0 < t) :
    ‖riemannZeta ((t : ℂ) * I)‖ ^ 2 * (2 * Real.pi * Real.cosh (Real.pi * t / 2))
      = t * Real.sinh (Real.pi * t / 2) * ‖riemannZeta (1 - (t : ℂ) * I)‖ ^ 2 := by
  have him : (1 - (t : ℂ) * I).im = -t := by simp
  have hs : ∀ n : ℕ, (1 - (t : ℂ) * I) ≠ -(n : ℂ) := by
    intro n h
    have := congrArg Complex.im h
    rw [him] at this; simp at this; linarith
  have hs' : (1 - (t : ℂ) * I) ≠ 1 := by
    intro h
    have := congrArg Complex.im h
    rw [him] at this; simp at this; linarith
  have fe := riemannZeta_one_sub hs hs'
  rw [show (1 : ℂ) - (1 - (t : ℂ) * I) = (t : ℂ) * I by ring] at fe
  have hcos : Complex.cos (↑Real.pi * (1 - (t : ℂ) * I) / 2)
      = ((Real.sinh (Real.pi * t / 2) : ℝ) : ℂ) * I := by
    rw [show (↑Real.pi : ℂ) * (1 - (t : ℂ) * I) / 2
        = ↑Real.pi / 2 - ((Real.pi * t / 2 : ℝ) : ℂ) * I by push_cast; ring,
      Complex.cos_pi_div_two_sub, Complex.sin_mul_I, Complex.ofReal_sinh]
  have h2π : (0 : ℝ) < 2 * Real.pi := by positivity
  have hpow : ‖(2 * (↑Real.pi : ℂ)) ^ (-(1 - (t : ℂ) * I))‖ = (2 * Real.pi)⁻¹ := by
    rw [show (2 * (↑Real.pi : ℂ)) = ((2 * Real.pi : ℝ) : ℂ) by push_cast; ring,
      Complex.norm_cpow_eq_rpow_re_of_pos h2π]
    have : (-(1 - (t : ℂ) * I)).re = -1 := by simp
    rw [this, Real.rpow_neg_one]
  have hsh : 0 < Real.sinh (Real.pi * t / 2) := Real.sinh_pos_iff.mpr (by positivity)
  have hch : 0 < Real.cosh (Real.pi * t / 2) := Real.cosh_pos _
  have hn : ‖riemannZeta ((t : ℂ) * I)‖
      = 2 * (2 * Real.pi)⁻¹ * ‖Complex.Gamma (1 - (t : ℂ) * I)‖ * Real.sinh (Real.pi * t / 2)
        * ‖riemannZeta (1 - (t : ℂ) * I)‖ := by
    rw [fe, norm_mul, norm_mul, norm_mul, norm_mul, hpow, hcos, norm_mul, Complex.norm_real,
      Complex.norm_I, mul_one, Real.norm_of_nonneg hsh.le]
    norm_num
  have hG := norm_Gamma_one_sub_mul_I_sq t ht
  have hsinh2 : Real.sinh (Real.pi * t)
      = 2 * Real.sinh (Real.pi * t / 2) * Real.cosh (Real.pi * t / 2) := by
    rw [← Real.sinh_two_mul]; ring_nf
  rw [hsinh2] at hG
  rw [hn]
  set G := ‖Complex.Gamma (1 - (t : ℂ) * I)‖
  set Z := ‖riemannZeta (1 - (t : ℂ) * I)‖
  set s := Real.sinh (Real.pi * t / 2)
  set c := Real.cosh (Real.pi * t / 2)
  have hpi : Real.pi ≠ 0 := Real.pi_ne_zero
  have expand : (2 * (2 * Real.pi)⁻¹ * G * s * Z) ^ 2 * (2 * Real.pi * c)
      = G ^ 2 * s ^ 2 * Z ^ 2 * (2 * c / Real.pi) := by
    field_simp
  rw [expand, hG]
  field_simp

theorem norm_zeta_mul_I_sq_le (t : ℝ) (ht : 0 < t) :
    ‖riemannZeta ((t : ℂ) * I)‖ ^ 2 ≤ t / (2 * Real.pi) * ‖riemannZeta (1 - (t : ℂ) * I)‖ ^ 2 := by
  have h := norm_zeta_mul_I_sq t ht
  have hch : 0 < Real.cosh (Real.pi * t / 2) := Real.cosh_pos _
  have hsc : Real.sinh (Real.pi * t / 2) ≤ Real.cosh (Real.pi * t / 2) := by
    rw [Real.sinh_eq, Real.cosh_eq]
    have := Real.exp_pos (-(Real.pi * t / 2)); linarith
  have hden : 0 < 2 * Real.pi * Real.cosh (Real.pi * t / 2) := by positivity
  rw [div_mul_eq_mul_div, le_div_iff₀ (by positivity : (0 : ℝ) < 2 * Real.pi)]
  have hZ : 0 ≤ ‖riemannZeta (1 - (t : ℂ) * I)‖ ^ 2 := by positivity
  have hzeta : 0 ≤ ‖riemannZeta ((t : ℂ) * I)‖ ^ 2 := by positivity
  nlinarith [mul_le_mul_of_nonneg_left hsc (mul_nonneg ht.le hZ)]

end ZetaImagAxis

#print axioms ZetaImagAxis.norm_Gamma_one_sub_mul_I_sq
#print axioms ZetaImagAxis.norm_zeta_mul_I_sq
#print axioms ZetaImagAxis.norm_zeta_mul_I_sq_le
