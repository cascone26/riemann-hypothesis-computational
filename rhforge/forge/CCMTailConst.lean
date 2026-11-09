import CCMCauchyLowner

/-!
# S2974 — far zeros contribute an (almost) CONSTANT rank-one block to the CCM zero side (lane 2)

Numerics (Report, HANDLE-grade cross-check, python parts(λ²=20,28; N=12) + first K zeros, w = ccmW = 4sin²(γL/2)/L,
C_{nγ} = 1/(2πn/L − γ), γ over ±zeros): A − G_K ≈ μ_K·𝟙𝟙ᵀ with entries ≈ 6.1e-3 (K=50), 3.3e-3 (K=150), spread 1e-5,
and the Löwner generator of A equals π·C·w up to a slope matching the Riemann–von Mangoldt tail to <1%
(K=30: −0.01233 vs −0.01188 … K=150: −0.004908 vs −0.004918).  So the forge zero-side model reproduces A
(the explicit formula in Löwner form).  This brick proves WHY the tail is ≈ constant:
* `cauchy_pair_near_const` — |a|,|b| ≤ V, 2V ≤ |x| ⇒ |1/((a−x)(b−x)) − 1/x²| ≤ 10V/|x|³.
* `cauchy_tail_entry_near_const` — for a weighted far-zero block (all |x_j| ≥ 2V, w ≥ 0):
  |Σ_j w_j/((ν_i−x_j)(ν_k−x_j)) − Σ_j w_j/x_j²| ≤ 10V·Σ_j w_j/|x_j|³.
-/

open Finset

namespace CCMTail

theorem cauchy_pair_near_const {a b x V : ℝ} (ha : |a| ≤ V) (hb : |b| ≤ V) (hx : 2 * V ≤ |x|)
    (hx0 : 0 < |x|) :
    |1 / ((a - x) * (b - x)) - 1 / x ^ 2| ≤ 10 * V / |x| ^ 3 := by
  have hV : 0 ≤ V := le_trans (abs_nonneg a) ha
  have hax : |x| / 2 ≤ |a - x| := by
    have := abs_sub_abs_le_abs_sub x a
    rw [abs_sub_comm] at this; linarith
  have hbx : |x| / 2 ≤ |b - x| := by
    have := abs_sub_abs_le_abs_sub x b
    rw [abs_sub_comm] at this; linarith
  have hapos : 0 < |a - x| := by linarith
  have hbpos : 0 < |b - x| := by linarith
  have hane : a - x ≠ 0 := abs_pos.mp hapos
  have hbne : b - x ≠ 0 := abs_pos.mp hbpos
  have hxne : x ≠ 0 := abs_pos.mp hx0
  have heq : 1 / ((a - x) * (b - x)) - 1 / x ^ 2
      = (x * (a + b) - a * b) / (x ^ 2 * ((a - x) * (b - x))) := by
    field_simp; ring
  rw [heq, abs_div]
  have hnum : |x * (a + b) - a * b| ≤ 5 * V * |x| / 2 := by
    calc |x * (a + b) - a * b| ≤ |x * (a + b)| + |a * b| := abs_sub _ _
      _ = |x| * |a + b| + |a| * |b| := by rw [abs_mul, abs_mul]
      _ ≤ |x| * (V + V) + V * V := by
          gcongr
          exact le_trans (abs_add_le a b) (add_le_add ha hb)
      _ ≤ 5 * V * |x| / 2 := by nlinarith
  have hden : |x| ^ 4 / 4 ≤ |x ^ 2 * ((a - x) * (b - x))| := by
    rw [abs_mul, abs_mul, abs_pow]
    have : |x| / 2 * (|x| / 2) ≤ |a - x| * |b - x| :=
      mul_le_mul hax hbx (by positivity) hapos.le
    nlinarith [sq_nonneg |x|]
  have hdpos : 0 < |x| ^ 4 / 4 := by positivity
  calc |x * (a + b) - a * b| / |x ^ 2 * ((a - x) * (b - x))|
      ≤ (5 * V * |x| / 2) / (|x| ^ 4 / 4) := by
        apply div_le_div₀ (by positivity) hnum hdpos hden
    _ = 10 * V / |x| ^ 3 := by field_simp; ring

theorem cauchy_tail_entry_near_const {m : Type*} [Fintype m] {V : ℝ} (a b : ℝ) (x w : m → ℝ)
    (ha : |a| ≤ V) (hb : |b| ≤ V) (hx : ∀ j, 2 * V ≤ |x j|) (hx0 : ∀ j, 0 < |x j|)
    (hw : ∀ j, 0 ≤ w j) :
    |∑ j, w j * (1 / ((a - x j) * (b - x j))) - ∑ j, w j * (1 / x j ^ 2)|
      ≤ 10 * V * ∑ j, w j / |x j| ^ 3 := by
  rw [← Finset.sum_sub_distrib, Finset.mul_sum]
  refine le_trans (Finset.abs_sum_le_sum_abs _ _) (Finset.sum_le_sum (fun j _ => ?_))
  rw [← mul_sub, abs_mul, abs_of_nonneg (hw j)]
  have h := cauchy_pair_near_const ha hb (hx j) (hx0 j)
  calc w j * |1 / ((a - x j) * (b - x j)) - 1 / x j ^ 2|
      ≤ w j * (10 * V / |x j| ^ 3) := mul_le_mul_of_nonneg_left h (hw j)
    _ = 10 * V * (w j / |x j| ^ 3) := by ring

end CCMTail

#print axioms CCMTail.cauchy_pair_near_const
#print axioms CCMTail.cauchy_tail_entry_near_const
