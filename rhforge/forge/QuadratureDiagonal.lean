/-
QuadratureDiagonal.lean (HP RH center, LEDGER S2896, 2026-11-07)
Generic linear-algebra form of S2868 (the Gauss-quadrature face of CCM 2511.22755 Lemma 5.4):
for a rank-2 displacement-structured symmetric T written as T = C A Cᵀ over Cauchy-type columns
(DC + CΘ = J), if some u with Cᵀu = 0 has nonzero sum (Shi's `he_sum` / CCM's ⟨δ|ξ⟩ = 1), then
A commutes with Θ — hence A is diagonal when Θ has distinct diagonal entries.
-/
import Mathlib

open Matrix

namespace QuadDiag

variable {m k : Type*} [Fintype m] [Fintype k] [DecidableEq m] [DecidableEq k]

lemma vecMulVec_mul' {p q r : Type*} [Fintype q] (x : p → ℝ) (y : q → ℝ) (M : Matrix q r ℝ) :
    vecMulVec x y * M = vecMulVec x (Mᵀ *ᵥ y) := by
  ext i j; simp [Matrix.mul_apply, vecMulVec_apply, mulVec, dotProduct, Finset.mul_sum, mul_comm, mul_left_comm]

lemma mul_vecMulVec' {p q r : Type*} [Fintype q] (M : Matrix p q ℝ) (x : q → ℝ) (y : r → ℝ) :
    M * vecMulVec x y = vecMulVec (M *ᵥ x) y := by
  ext i j; simp [Matrix.mul_apply, vecMulVec_apply, mulVec, dotProduct, Finset.sum_mul, mul_assoc]

lemma vecMulVec_mulVec' {p q : Type*} [Fintype q] (x : p → ℝ) (y u : q → ℝ) :
    vecMulVec x y *ᵥ u = (y ⬝ᵥ u) • x := by
  ext i; simp [mulVec, vecMulVec_apply, dotProduct, Finset.mul_sum, mul_comm, mul_left_comm, mul_assoc]

/-- Step 3: C X Cᵀ = 0 with C injective forces X = 0. -/
lemma eq_zero_of_sandwich (C : Matrix m k ℝ) (X : Matrix k k ℝ)
    (hinj : Function.Injective C.mulVec) (h : C * X * Cᵀ = 0) : X = 0 := by
  have h1 : ∀ v, X *ᵥ (Cᵀ *ᵥ v) = 0 := by
    intro v; apply hinj
    rw [Matrix.mulVec_zero, Matrix.mulVec_mulVec, Matrix.mulVec_mulVec, h, zero_mulVec]
  have h2 : C * Xᵀ = 0 := by
    have hXC : X * Cᵀ = 0 := by
      ext a i
      have := congrFun (h1 (Pi.single i 1)) a
      rw [Matrix.mulVec_mulVec, mulVec_single_one] at this
      simpa [transpose_apply] using this
    have := congrArg transpose hXC; simpa [transpose_mul] using this
  ext a b
  have h3 : Xᵀ *ᵥ Pi.single a 1 = 0 := by
    apply hinj; rw [Matrix.mulVec_zero, Matrix.mulVec_mulVec, h2, zero_mulVec]
  have := congrFun h3 b
  rw [mulVec_single_one] at this
  simpa [transpose_apply] using this

/-- **Quadrature commutation theorem.** -/
theorem theta_comm (T D : Matrix m m ℝ) (Θ A : Matrix k k ℝ) (C : Matrix m k ℝ) (b : m → ℝ)
    (hD : Dᵀ = D) (hΘ : Θᵀ = Θ) (hA : Aᵀ = A)
    (hdisp : D * T - T * D = vecMulVec b (fun _ => 1) - vecMulVec (fun _ => 1) b)
    (hT : T = C * A * Cᵀ)
    (hC : D * C + C * Θ = vecMulVec (fun _ => (1:ℝ)) (fun _ => 1))
    (hinj : Function.Injective C.mulVec)
    (hu : ∃ u : m → ℝ, Cᵀ *ᵥ u = 0 ∧ ∑ i, u i ≠ 0) :
    A * Θ = Θ * A := by
  set J : Matrix m k ℝ := vecMulVec (fun _ => (1:ℝ)) (fun _ => 1) with hJ
  have hDC : D * C = J - C * Θ := by rw [← hC]; abel
  have hCD : Cᵀ * D = Jᵀ - Θ * Cᵀ := by
    have := congrArg transpose hDC
    rw [transpose_mul, hD, transpose_sub, transpose_mul, hΘ] at this
    exact this
  set y : m → ℝ := (C * A) *ᵥ (fun _ => 1) with hy
  have hJA : J * A * Cᵀ = vecMulVec (fun _ => 1) y := by
    rw [hJ, vecMulVec_mul', vecMulVec_mul', hy, transpose_transpose, hA, Matrix.mulVec_mulVec]
  have hAJ : C * A * Jᵀ = vecMulVec y (fun _ => 1) := by
    rw [hJ, transpose_vecMulVec, mul_vecMulVec', hy]
  set X : Matrix k k ℝ := A * Θ - Θ * A
  have hkey : C * X * Cᵀ = vecMulVec (b + y) (fun _ => 1) - vecMulVec (fun _ => 1) (b + y) := by
    have e1 : D * T = J * A * Cᵀ - C * Θ * A * Cᵀ := by
      rw [hT, ← Matrix.mul_assoc, ← Matrix.mul_assoc, hDC, Matrix.sub_mul, Matrix.sub_mul]
    have e2 : T * D = C * A * Jᵀ - C * A * Θ * Cᵀ := by
      rw [hT, Matrix.mul_assoc, hCD, Matrix.mul_sub]; simp [Matrix.mul_assoc]
    have e3 : C * X * Cᵀ = (D * T - T * D) - J * A * Cᵀ + C * A * Jᵀ := by
      rw [e1, e2]; simp only [X, Matrix.mul_sub, Matrix.sub_mul, Matrix.mul_assoc]; abel
    rw [e3, hdisp, hJA, hAJ]
    ext i j; simp [vecMulVec_apply, Pi.add_apply]; ring
  obtain ⟨u, hCu, hsum⟩ := hu
  set w := b + y
  have hw : (∑ i, u i) • w = (w ⬝ᵥ u) • (fun _ => (1:ℝ)) := by
    have h0 : (C * X * Cᵀ) *ᵥ u = (vecMulVec w (fun _ => 1) - vecMulVec (fun _ => 1) w) *ᵥ u := by
      rw [hkey]
    have lhs0 : (C * X * Cᵀ) *ᵥ u = 0 := by
      rw [← Matrix.mulVec_mulVec, hCu, Matrix.mulVec_zero]
    rw [lhs0, Matrix.sub_mulVec, vecMulVec_mulVec', vecMulVec_mulVec'] at h0
    have h1 : ((fun _ => (1:ℝ)) ⬝ᵥ u) = ∑ i, u i := by simp [dotProduct]
    rw [h1] at h0
    exact (sub_eq_zero.mp h0.symm)
  have hwc : w = ((w ⬝ᵥ u) / ∑ i, u i) • (fun _ => (1:ℝ)) := by
    rw [div_eq_inv_mul, mul_smul, ← hw, smul_smul, inv_mul_cancel₀ hsum, one_smul]
  have hzero : C * X * Cᵀ = 0 := by
    rw [hkey, hwc]; ext i j; simp [vecMulVec_apply]
  have := eq_zero_of_sandwich C X hinj hzero
  exact sub_eq_zero.mp this

/-- Diagonality: if Θ is diagonal with distinct entries, A is diagonal. -/
theorem A_diagonal (Θ A : Matrix k k ℝ) (θ : k → ℝ) (hΘd : Θ = diagonal θ) (hθ : Function.Injective θ)
    (hcomm : A * Θ = Θ * A) {i j : k} (hij : i ≠ j) : A i j = 0 := by
  have := congrFun (congrFun hcomm i) j
  rw [hΘd, mul_diagonal, diagonal_mul] at this
  have hne : θ j - θ i ≠ 0 := sub_ne_zero.mpr (fun h => hij (hθ h.symm))
  have : A i j * (θ j - θ i) = 0 := by linarith [this]
  rcases mul_eq_zero.mp this with h | h
  · exact h
  · exact absurd h hne

end QuadDiag

#print axioms QuadDiag.theta_comm
#print axioms QuadDiag.A_diagonal
