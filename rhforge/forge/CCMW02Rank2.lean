import Mathlib

/-!
# S2966 — arithmetic side, first brick: the CCM pole block W02 is rank 2 of signature (1,1) (lane 2)

From `operators/ccm_parts.py`: on frequencies ν : ι → ℝ (ν = n ∈ {−N..N}),
  `W02[i,j] = 32 L sinh²(L/4) (L² − 16π² ν_i ν_j) / ((L² + 16π² ν_i²)(L² + 16π² ν_j²))`.
With `u_i = L/(L²+16π²ν_i²)`, `v_i = 4πν_i/(L²+16π²ν_i²)`, `c = 32 L sinh²(L/4)`:
* `W02_eq` — `W02 = c • (u uᵀ − v vᵀ)`.
* `quad_vecMulVec_self` — `xᵀ(w wᵀ)x = (w·x)²`.
* `W02_quad` — `xᵀ W02 x = c((u·x)² − (v·x)²)`.
* `W02_quad_le` — `xᵀ W02 x ≤ c (u·x)²` (L>0); `W02_quad_nonpos_of_perp` — `u·x = 0 ⇒ xᵀ W02 x ≤ 0`.
Reading (Report): the pole block contributes positivity in exactly ONE direction (u); on u^⊥ it is ≤ 0, so
positivity of W02 − W_R − W_p on the codimension-1 space u^⊥ must come entirely from −(W_R + W_p).
Combined with S2964 (zero side shrinks at most linearly in the window), this localises the e^{−κλ²} margin
to the archimedean + prime blocks restricted to u^⊥ — the target of the next arithmetic-side bricks.
No novelty claimed (W02 is the pole term of the Weil form).
-/

open Matrix Finset Real

namespace CCMW02

variable {ι : Type*} [Fintype ι]

noncomputable def W02 (L : ℝ) (ν : ι → ℝ) : Matrix ι ι ℝ := fun i j =>
  32 * L * Real.sinh (L / 4) ^ 2 * (L ^ 2 - 16 * π ^ 2 * ν i * ν j)
    / ((L ^ 2 + 16 * π ^ 2 * ν i ^ 2) * (L ^ 2 + 16 * π ^ 2 * ν j ^ 2))

noncomputable def uW (L : ℝ) (ν : ι → ℝ) : ι → ℝ := fun i => L / (L ^ 2 + 16 * π ^ 2 * ν i ^ 2)
noncomputable def vW (L : ℝ) (ν : ι → ℝ) : ι → ℝ := fun i =>
  4 * π * ν i / (L ^ 2 + 16 * π ^ 2 * ν i ^ 2)
noncomputable def cW (L : ℝ) : ℝ := 32 * L * Real.sinh (L / 4) ^ 2

theorem den_pos {L : ℝ} (hL : 0 < L) (y : ℝ) : 0 < L ^ 2 + 16 * π ^ 2 * y ^ 2 := by positivity

omit [Fintype ι] in
theorem W02_eq {L : ℝ} (hL : 0 < L) (ν : ι → ℝ) :
    W02 L ν = cW L • (vecMulVec (uW L ν) (uW L ν) - vecMulVec (vW L ν) (vW L ν)) := by
  ext i j
  have hi := (den_pos hL (ν i)).ne'
  have hj := (den_pos hL (ν j)).ne'
  simp only [W02, uW, vW, cW, Matrix.smul_apply, Matrix.sub_apply, vecMulVec_apply, smul_eq_mul]
  field_simp
  ring

theorem quad_vecMulVec_self (w x : ι → ℝ) : x ⬝ᵥ (vecMulVec w w *ᵥ x) = (w ⬝ᵥ x) ^ 2 := by
  simp only [dotProduct, mulVec, vecMulVec_apply]
  rw [sq, Finset.sum_mul_sum]
  apply Finset.sum_congr rfl
  intro i _
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro j _
  ring

theorem W02_quad {L : ℝ} (hL : 0 < L) (ν x : ι → ℝ) :
    x ⬝ᵥ (W02 L ν *ᵥ x) = cW L * ((uW L ν ⬝ᵥ x) ^ 2 - (vW L ν ⬝ᵥ x) ^ 2) := by
  rw [W02_eq hL, Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul, Matrix.sub_mulVec,
    dotProduct_sub, quad_vecMulVec_self, quad_vecMulVec_self]

theorem cW_nonneg {L : ℝ} (hL : 0 < L) : 0 ≤ cW L := by unfold cW; positivity

theorem W02_quad_le {L : ℝ} (hL : 0 < L) (ν x : ι → ℝ) :
    x ⬝ᵥ (W02 L ν *ᵥ x) ≤ cW L * (uW L ν ⬝ᵥ x) ^ 2 := by
  rw [W02_quad hL]
  have := cW_nonneg hL
  nlinarith [sq_nonneg (vW L ν ⬝ᵥ x)]

theorem W02_quad_nonpos_of_perp {L : ℝ} (hL : 0 < L) (ν x : ι → ℝ) (hx : uW L ν ⬝ᵥ x = 0) :
    x ⬝ᵥ (W02 L ν *ᵥ x) ≤ 0 := by
  have h := W02_quad_le hL ν x
  rw [hx] at h
  simpa using h

theorem W02_quad_nonneg_of_perp_v {L : ℝ} (hL : 0 < L) (ν x : ι → ℝ) (hx : vW L ν ⬝ᵥ x = 0) :
    0 ≤ x ⬝ᵥ (W02 L ν *ᵥ x) := by
  rw [W02_quad hL, hx]
  have := cW_nonneg hL
  simp only [ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, sub_zero]
  positivity

end CCMW02

#print axioms CCMW02.W02_eq
#print axioms CCMW02.quad_vecMulVec_self
#print axioms CCMW02.W02_quad
#print axioms CCMW02.W02_quad_le
#print axioms CCMW02.W02_quad_nonpos_of_perp
#print axioms CCMW02.W02_quad_nonneg_of_perp_v
