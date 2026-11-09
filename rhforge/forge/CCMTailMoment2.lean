import CCMTailPerp

/-!
# S2983 — moment hierarchy, level 2: on span{𝟙, ν, ν²}^⊥ the paired far tail is O(V⁴Σw/x⁶) (lane 2)

S2980 (negative): the CCM margin lives in the far tail's fine structure.  Expanding the ±x pair kernel:
  `pair(a,b,x) = 2/x² + 2(ab + a² + b²)/x⁴ + R`,  `R = (2S(a²+b²)x² − 2a²b²x² − 2Sa²b²)/(x⁴(x²−a²)(x²−b²))`, S = ab+a²+b².
The first two orders are a rank-≤4 matrix in span{𝟙, ν, ν²}; on its orthogonal complement only R survives.
* `pair_second_order` — |a|,|b| ≤ V, 2V ≤ |x| ⇒ |pair − 2/x² − 2(ab+a²+b²)/x⁴| ≤ 29V⁴/x⁶.
* `quad_vecMulVec` — xᵀ(u vᵀ)x = (u·x)(v·x).
* `near_matrix_quad` — entries of M within δ of P ⇒ |xᵀMx − xᵀPx| ≤ δ(Σ|x|)².
* `moment2_quad` — P = α𝟙𝟙ᵀ + β(ννᵀ + ν²𝟙ᵀ + 𝟙ν²ᵀ), x ⊥ 𝟙, ν, ν² ⇒ |xᵀMx| ≤ δ(Σ|x|)².
-/

open Finset Matrix

namespace CCMMoment2

set_option maxHeartbeats 1000000 in
theorem pair_second_order {a b x V : ℝ} (ha : |a| ≤ V) (hb : |b| ≤ V) (hx : 2 * V ≤ |x|)
    (hx0 : x ≠ 0) :
    |1 / ((a - x) * (b - x)) + 1 / ((a + x) * (b + x)) - 2 / x ^ 2
        - 2 * (a * b + a ^ 2 + b ^ 2) / x ^ 4| ≤ 29 * V ^ 4 / x ^ 6 := by
  have hV : 0 ≤ V := le_trans (abs_nonneg a) ha
  have hA : a ^ 2 ≤ V ^ 2 := by rw [← sq_abs a]; exact pow_le_pow_left₀ (abs_nonneg a) ha 2
  have hB : b ^ 2 ≤ V ^ 2 := by rw [← sq_abs b]; exact pow_le_pow_left₀ (abs_nonneg b) hb 2
  have hX : 4 * V ^ 2 ≤ x ^ 2 := by rw [← sq_abs x]; nlinarith [abs_nonneg x]
  have hx2 : 0 < x ^ 2 := by positivity
  have hdA : 3 * x ^ 2 / 4 ≤ x ^ 2 - a ^ 2 := by linarith
  have hdB : 3 * x ^ 2 / 4 ≤ x ^ 2 - b ^ 2 := by linarith
  have hpA : 0 < x ^ 2 - a ^ 2 := by linarith
  have hpB : 0 < x ^ 2 - b ^ 2 := by linarith
  have h1 : a - x ≠ 0 := by
    intro h; have : a = x := by linarith
    rw [this] at hpA; linarith
  have h2 : b - x ≠ 0 := by
    intro h; have : b = x := by linarith
    rw [this] at hpB; linarith
  have h3 : a + x ≠ 0 := by
    intro h; have : a = -x := by linarith
    rw [this] at hpA; nlinarith
  have h4 : b + x ≠ 0 := by
    intro h; have : b = -x := by linarith
    rw [this] at hpB; nlinarith
  set S := a * b + a ^ 2 + b ^ 2 with hSdef
  have heq : 1 / ((a - x) * (b - x)) + 1 / ((a + x) * (b + x)) - 2 / x ^ 2 - 2 * S / x ^ 4
      = (2 * S * (a ^ 2 + b ^ 2) * x ^ 2 - 2 * a ^ 2 * b ^ 2 * x ^ 2 - 2 * S * a ^ 2 * b ^ 2)
        / (x ^ 4 * ((x ^ 2 - a ^ 2) * (x ^ 2 - b ^ 2))) := by
    rw [hSdef]; field_simp; ring
  rw [heq]
  have hDpos : 0 < x ^ 4 * ((x ^ 2 - a ^ 2) * (x ^ 2 - b ^ 2)) := by positivity
  rw [abs_div, abs_of_pos hDpos]
  have hS0 : 0 ≤ S := by rw [hSdef]; nlinarith [sq_nonneg (a + b)]
  have hS : S ≤ 3 * V ^ 2 := by rw [hSdef]; nlinarith [sq_nonneg (a - b)]
  have hab : a ^ 2 * b ^ 2 ≤ V ^ 2 * V ^ 2 := mul_le_mul hA hB (sq_nonneg b) (sq_nonneg V)
  have hab0 : 0 ≤ a ^ 2 * b ^ 2 := by positivity
  have hsum : a ^ 2 + b ^ 2 ≤ 2 * V ^ 2 := by linarith
  have hV2 : V ^ 2 ≤ x ^ 2 / 4 := by linarith
  have t1 : 0 ≤ 2 * S * (a ^ 2 + b ^ 2) * x ^ 2 := by positivity
  have t1' : 2 * S * (a ^ 2 + b ^ 2) * x ^ 2 ≤ 12 * V ^ 4 * x ^ 2 := by
    have : S * (a ^ 2 + b ^ 2) ≤ 3 * V ^ 2 * (2 * V ^ 2) :=
      mul_le_mul hS hsum (by positivity) (by positivity)
    nlinarith
  have t2 : 2 * a ^ 2 * b ^ 2 * x ^ 2 ≤ 2 * V ^ 4 * x ^ 2 := by nlinarith
  have t3 : 2 * S * a ^ 2 * b ^ 2 ≤ 6 * V ^ 4 * (x ^ 2 / 4) := by
    have m1 : S * (a ^ 2 * b ^ 2) ≤ 3 * V ^ 2 * (V ^ 2 * V ^ 2) :=
      mul_le_mul hS hab hab0 (by positivity)
    have m2 : V ^ 4 * V ^ 2 ≤ V ^ 4 * (x ^ 2 / 4) :=
      mul_le_mul_of_nonneg_left hV2 (by positivity)
    have e1 : 3 * V ^ 2 * (V ^ 2 * V ^ 2) = 3 * (V ^ 4 * V ^ 2) := by ring
    have e2 : 2 * S * a ^ 2 * b ^ 2 = 2 * (S * (a ^ 2 * b ^ 2)) := by ring
    rw [e2]; linarith
  have t3' : 0 ≤ 2 * S * a ^ 2 * b ^ 2 := by positivity
  have hnum : |2 * S * (a ^ 2 + b ^ 2) * x ^ 2 - 2 * a ^ 2 * b ^ 2 * x ^ 2 - 2 * S * a ^ 2 * b ^ 2|
      ≤ 16 * V ^ 4 * x ^ 2 := by
    have t2' : 0 ≤ 2 * a ^ 2 * b ^ 2 * x ^ 2 := by positivity
    have hx2' : 0 ≤ V ^ 4 * x ^ 2 := by positivity
    rw [abs_le]; constructor <;> linarith
  have hden : 9 * x ^ 8 / 16 ≤ x ^ 4 * ((x ^ 2 - a ^ 2) * (x ^ 2 - b ^ 2)) := by
    have hprod : (3 * x ^ 2 / 4) * (3 * x ^ 2 / 4) ≤ (x ^ 2 - a ^ 2) * (x ^ 2 - b ^ 2) :=
      mul_le_mul hdA hdB (by positivity) hpA.le
    have hx4 : 0 ≤ x ^ 4 := by positivity
    nlinarith
  have hd0 : 0 < 9 * x ^ 8 / 16 := by positivity
  calc |2 * S * (a ^ 2 + b ^ 2) * x ^ 2 - 2 * a ^ 2 * b ^ 2 * x ^ 2 - 2 * S * a ^ 2 * b ^ 2|
        / (x ^ 4 * ((x ^ 2 - a ^ 2) * (x ^ 2 - b ^ 2)))
      ≤ (16 * V ^ 4 * x ^ 2) / (9 * x ^ 8 / 16) := div_le_div₀ (by positivity) hnum hd0 hden
    _ = 256 / 9 * V ^ 4 / x ^ 6 := by field_simp; ring
    _ ≤ 29 * V ^ 4 / x ^ 6 := by
        apply div_le_div_of_nonneg_right _ (by positivity)
        nlinarith [pow_nonneg hV 4]

variable {n : Type*} [Fintype n]

theorem quad_vecMulVec (u v x : n → ℝ) : x ⬝ᵥ (vecMulVec u v *ᵥ x) = (u ⬝ᵥ x) * (v ⬝ᵥ x) := by
  simp only [dotProduct, mulVec, vecMulVec_apply]
  rw [Finset.sum_mul_sum]
  apply Finset.sum_congr rfl; intro i _
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl; intro j _
  ring

theorem near_matrix_quad (M P : Matrix n n ℝ) (δ : ℝ) (h : ∀ i k, |M i k - P i k| ≤ δ)
    (x : n → ℝ) : |x ⬝ᵥ (M *ᵥ x) - x ⬝ᵥ (P *ᵥ x)| ≤ δ * (∑ i, |x i|) ^ 2 := by
  have hq := CCMTailQuad.near_const_quad (M - P) 0 δ (fun i k => by simpa using h i k) x
  rw [Matrix.sub_mulVec, dotProduct_sub] at hq
  simpa using hq

theorem moment2_quad (M : Matrix n n ℝ) (ν : n → ℝ) (α β δ : ℝ)
    (h : ∀ i k, |M i k - (α + β * (ν i * ν k + ν i ^ 2 + ν k ^ 2))| ≤ δ) (x : n → ℝ)
    (h1 : (fun _ => (1 : ℝ)) ⬝ᵥ x = 0) (hν : ν ⬝ᵥ x = 0) (hν2 : (fun i => ν i ^ 2) ⬝ᵥ x = 0) :
    |x ⬝ᵥ (M *ᵥ x)| ≤ δ * (∑ i, |x i|) ^ 2 := by
  set P : Matrix n n ℝ := α • vecMulVec (fun _ => (1 : ℝ)) (fun _ => 1)
    + β • (vecMulVec ν ν + vecMulVec (fun i => ν i ^ 2) (fun _ => 1)
      + vecMulVec (fun _ => (1 : ℝ)) (fun i => ν i ^ 2)) with hP
  have hPe : ∀ i k, P i k = α + β * (ν i * ν k + ν i ^ 2 + ν k ^ 2) := by
    intro i k; simp [hP, vecMulVec_apply]; ring
  have hq := near_matrix_quad M P δ (fun i k => by rw [hPe]; exact h i k) x
  have hP0 : x ⬝ᵥ (P *ᵥ x) = 0 := by
    simp only [hP, Matrix.add_mulVec, Matrix.smul_mulVec, dotProduct_add, dotProduct_smul,
      smul_eq_mul, quad_vecMulVec, h1, hν, hν2]
    ring
  rw [hP0, sub_zero] at hq
  exact hq

end CCMMoment2

#print axioms CCMMoment2.pair_second_order
#print axioms CCMMoment2.quad_vecMulVec
#print axioms CCMMoment2.near_matrix_quad
#print axioms CCMMoment2.moment2_quad
