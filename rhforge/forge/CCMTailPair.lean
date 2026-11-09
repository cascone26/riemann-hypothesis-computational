import CCMTailQuad

/-!
# S2978 — ±γ pairing kills the odd term: far-zero tail error O(V²/γ⁴) (lane 2)

S2976's tail bound used `|1/((a−x)(b−x)) − 1/x²| ≤ 10V/|x|³` and was ~60× loose (δ > μ at λ²=20).  Zeros come in pairs
±γ with equal weight (ccmW is even in γ), and for a pair the (a+b)/x³ term cancels:
  `1/((a−x)(b−x)) + 1/((a+x)(b+x)) − 2/x² = (2x²(ab+a²+b²) − 2a²b²) / (x²(x²−a²)(x²−b²))`.
* `pair_near_const` — |a|,|b| ≤ V, 2V ≤ |x| ⇒ |pair − 2/x²| ≤ 11V²/x⁴.
* `pair_tail_entry` — Σ_j w_j·pair_j within 11V²·Σ w_j/x_j⁴ of Σ 2w_j/x_j².
* `pair_tail_quad` — |xᵀT x − μ(Σx)²| ≤ δ(Σ|x|)² with μ = Σ 2w_j/x_j², δ = 11V²Σ w_j/x_j⁴.
-/

open Finset Matrix

namespace CCMTailPair

theorem pair_near_const {a b x V : ℝ} (ha : |a| ≤ V) (hb : |b| ≤ V) (hx : 2 * V ≤ |x|) (hx0 : x ≠ 0) :
    |1 / ((a - x) * (b - x)) + 1 / ((a + x) * (b + x)) - 2 / x ^ 2| ≤ 11 * V ^ 2 / x ^ 4 := by
  have hV : 0 ≤ V := le_trans (abs_nonneg a) ha
  have hA : a ^ 2 ≤ V ^ 2 := by rw [← sq_abs a]; exact pow_le_pow_left₀ (abs_nonneg a) ha 2
  have hB : b ^ 2 ≤ V ^ 2 := by rw [← sq_abs b]; exact pow_le_pow_left₀ (abs_nonneg b) hb 2
  have hX : 4 * V ^ 2 ≤ x ^ 2 := by
    rw [← sq_abs x]; nlinarith [abs_nonneg x]
  have hx2 : 0 < x ^ 2 := by positivity
  have hdA : 3 * x ^ 2 / 4 ≤ x ^ 2 - a ^ 2 := by linarith
  have hdB : 3 * x ^ 2 / 4 ≤ x ^ 2 - b ^ 2 := by linarith
  have hpA : 0 < x ^ 2 - a ^ 2 := by linarith
  have hpB : 0 < x ^ 2 - b ^ 2 := by linarith
  have h1 : a - x ≠ 0 := by intro h; have : a = x := by linarith
                            rw [this] at hpA; linarith
  have h2 : b - x ≠ 0 := by intro h; have : b = x := by linarith
                            rw [this] at hpB; linarith
  have h3 : a + x ≠ 0 := by intro h; have : a = -x := by linarith
                            rw [this] at hpA; nlinarith
  have h4 : b + x ≠ 0 := by intro h; have : b = -x := by linarith
                            rw [this] at hpB; nlinarith
  have heq : 1 / ((a - x) * (b - x)) + 1 / ((a + x) * (b + x)) - 2 / x ^ 2
      = (2 * x ^ 2 * (a * b + a ^ 2 + b ^ 2) - 2 * a ^ 2 * b ^ 2)
        / (x ^ 2 * ((x ^ 2 - a ^ 2) * (x ^ 2 - b ^ 2))) := by
    field_simp
    ring
  rw [heq]
  have hDpos : 0 < x ^ 2 * ((x ^ 2 - a ^ 2) * (x ^ 2 - b ^ 2)) := by positivity
  rw [abs_div, abs_of_pos hDpos]
  have hS0 : 0 ≤ a * b + a ^ 2 + b ^ 2 := by nlinarith [sq_nonneg (a + b)]
  have hS : a * b + a ^ 2 + b ^ 2 ≤ 3 * V ^ 2 := by nlinarith [sq_nonneg (a - b)]
  have hab : a ^ 2 * b ^ 2 ≤ V ^ 2 * V ^ 2 := mul_le_mul hA hB (sq_nonneg b) (sq_nonneg V)
  have hnum : |2 * x ^ 2 * (a * b + a ^ 2 + b ^ 2) - 2 * a ^ 2 * b ^ 2| ≤ 6 * V ^ 2 * x ^ 2 := by
    rw [abs_le]; constructor
    · nlinarith [sq_nonneg a, sq_nonneg b, mul_nonneg (sq_nonneg a) (sq_nonneg b)]
    · nlinarith [mul_nonneg (sq_nonneg a) (sq_nonneg b)]
  have hden : 9 * x ^ 6 / 16 ≤ x ^ 2 * ((x ^ 2 - a ^ 2) * (x ^ 2 - b ^ 2)) := by
    have hprod : (3 * x ^ 2 / 4) * (3 * x ^ 2 / 4) ≤ (x ^ 2 - a ^ 2) * (x ^ 2 - b ^ 2) :=
      mul_le_mul hdA hdB (by positivity) hpA.le
    nlinarith
  have hd0 : 0 < 9 * x ^ 6 / 16 := by positivity
  calc |2 * x ^ 2 * (a * b + a ^ 2 + b ^ 2) - 2 * a ^ 2 * b ^ 2|
        / (x ^ 2 * ((x ^ 2 - a ^ 2) * (x ^ 2 - b ^ 2)))
      ≤ (6 * V ^ 2 * x ^ 2) / (9 * x ^ 6 / 16) := div_le_div₀ (by positivity) hnum hd0 hden
    _ = 32 / 3 * V ^ 2 / x ^ 4 := by field_simp; ring
    _ ≤ 11 * V ^ 2 / x ^ 4 := by
        apply div_le_div_of_nonneg_right _ (by positivity)
        nlinarith [sq_nonneg V]

theorem pair_tail_entry {m : Type*} [Fintype m] {V : ℝ} (a b : ℝ) (x w : m → ℝ)
    (ha : |a| ≤ V) (hb : |b| ≤ V) (hx : ∀ j, 2 * V ≤ |x j|) (hx0 : ∀ j, x j ≠ 0)
    (hw : ∀ j, 0 ≤ w j) :
    |∑ j, w j * (1 / ((a - x j) * (b - x j)) + 1 / ((a + x j) * (b + x j)))
        - ∑ j, w j * (2 / x j ^ 2)|
      ≤ 11 * V ^ 2 * ∑ j, w j / x j ^ 4 := by
  rw [← Finset.sum_sub_distrib, Finset.mul_sum]
  refine le_trans (Finset.abs_sum_le_sum_abs _ _) (Finset.sum_le_sum (fun j _ => ?_))
  rw [← mul_sub, abs_mul, abs_of_nonneg (hw j)]
  have h := pair_near_const ha hb (hx j) (hx0 j)
  calc w j * |1 / ((a - x j) * (b - x j)) + 1 / ((a + x j) * (b + x j)) - 2 / x j ^ 2|
      ≤ w j * (11 * V ^ 2 / x j ^ 4) := mul_le_mul_of_nonneg_left h (hw j)
    _ = 11 * V ^ 2 * (w j / x j ^ 4) := by ring

theorem pair_tail_quad {n m : Type*} [Fintype n] [Fintype m] {V : ℝ} (ν : n → ℝ) (xz w : m → ℝ)
    (hν : ∀ i, |ν i| ≤ V) (hx : ∀ j, 2 * V ≤ |xz j|) (hx0 : ∀ j, xz j ≠ 0) (hw : ∀ j, 0 ≤ w j)
    (x : n → ℝ) :
    |x ⬝ᵥ ((fun i k => ∑ j, w j * (1 / ((ν i - xz j) * (ν k - xz j))
          + 1 / ((ν i + xz j) * (ν k + xz j)))) *ᵥ x)
        - (∑ j, w j * (2 / xz j ^ 2)) * (∑ i, x i) ^ 2|
      ≤ (11 * V ^ 2 * ∑ j, w j / xz j ^ 4) * (∑ i, |x i|) ^ 2 :=
  CCMTailQuad.near_const_quad _ _ _ (fun i k =>
    pair_tail_entry (ν i) (ν k) xz w (hν i) (hν k) hx hx0 hw) x

end CCMTailPair

#print axioms CCMTailPair.pair_near_const
#print axioms CCMTailPair.pair_tail_entry
#print axioms CCMTailPair.pair_tail_quad
