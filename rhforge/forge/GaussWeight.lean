import Mathlib

/-!
# S2986 — Gaussian-weighted polynomial suprema (lane 1: three-lines edge inputs a, b)

S2984 needs `a ≥ sup_t ‖E(it)‖e^{−(t−T)²}` and `b ≥ sup_t ‖E(1+it)‖e^{1−(t−T)²}` with the SAME exponent in T as the
edge bounds (losing exponent kills convexity).  Exact-exponent Gaussian moment bound:
* `one_add_abs_le_mul` — 1 + |t| ≤ (1 + |T|)(1 + |t − T|).
* `one_add_rpow_le_exp` — u ≥ 0, p ≥ 0 ⇒ (1 + u)^p ≤ e^{pu}.
* `gauss_poly_weight` — p ≥ 0 ⇒ (1 + |t|)^p · e^{−(t−T)²} ≤ (1 + |T|)^p · e^{p²/4}   (for ALL real t, T).
So a polynomial edge bound C(1+|t|)^p localises to C e^{p²/4}(1+|T|)^p — no exponent loss.
-/

namespace GaussWeight

theorem one_add_abs_le_mul (t T : ℝ) : 1 + |t| ≤ (1 + |T|) * (1 + |t - T|) := by
  have h1 : |t| ≤ |T| + |t - T| := by
    have := abs_add_le T (t - T); simpa using this
  have h2 : 0 ≤ |T| * |t - T| := mul_nonneg (abs_nonneg _) (abs_nonneg _)
  nlinarith

theorem one_add_rpow_le_exp {u p : ℝ} (hu : 0 ≤ u) (hp : 0 ≤ p) :
    (1 + u) ^ p ≤ Real.exp (p * u) := by
  have h1 : 1 + u ≤ Real.exp u := by linarith [Real.add_one_le_exp u]
  calc (1 + u) ^ p ≤ (Real.exp u) ^ p := Real.rpow_le_rpow (by linarith) h1 hp
    _ = Real.exp (p * u) := by rw [← Real.exp_mul]; ring_nf

theorem gauss_poly_weight (t T : ℝ) {p : ℝ} (hp : 0 ≤ p) :
    (1 + |t|) ^ p * Real.exp (-(t - T) ^ 2) ≤ (1 + |T|) ^ p * Real.exp (p ^ 2 / 4) := by
  set u := |t - T| with hu
  have hu0 : 0 ≤ u := abs_nonneg _
  have hsq : (t - T) ^ 2 = u ^ 2 := by rw [hu, sq_abs]
  have hA : (1 + |t|) ^ p ≤ ((1 + |T|) * (1 + u)) ^ p :=
    Real.rpow_le_rpow (by positivity) (one_add_abs_le_mul t T) hp
  rw [Real.mul_rpow (by positivity) (by positivity)] at hA
  have hB := one_add_rpow_le_exp hu0 hp
  have hE : Real.exp (p * u) * Real.exp (-(t - T) ^ 2) ≤ Real.exp (p ^ 2 / 4) := by
    rw [← Real.exp_add, hsq]
    apply Real.exp_le_exp.mpr
    nlinarith [sq_nonneg (u - p / 2)]
  have hT0 : 0 ≤ (1 + |T|) ^ p := by positivity
  have hg0 : 0 ≤ Real.exp (-(t - T) ^ 2) := (Real.exp_pos _).le
  calc (1 + |t|) ^ p * Real.exp (-(t - T) ^ 2)
      ≤ (1 + |T|) ^ p * (1 + u) ^ p * Real.exp (-(t - T) ^ 2) :=
        mul_le_mul_of_nonneg_right hA hg0
    _ ≤ (1 + |T|) ^ p * Real.exp (p * u) * Real.exp (-(t - T) ^ 2) := by
        gcongr
    _ = (1 + |T|) ^ p * (Real.exp (p * u) * Real.exp (-(t - T) ^ 2)) := by ring
    _ ≤ (1 + |T|) ^ p * Real.exp (p ^ 2 / 4) := mul_le_mul_of_nonneg_left hE hT0

end GaussWeight

#print axioms GaussWeight.one_add_abs_le_mul
#print axioms GaussWeight.one_add_rpow_le_exp
#print axioms GaussWeight.gauss_poly_weight
