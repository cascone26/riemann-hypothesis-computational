/-
CauchyDetStep.lean (HP RH center, LEDGER S2913, 2026-11-08) — Cauchy determinant, step 2 (one-step recursion).
det C_{n+1}(x,y) = (x_n+y_n)⁻¹ · Π_i (x_n−x_i)/(x_i+y_n) · det C_n(x',y') · Π_j (y_n−y_j)/(x_n+y_j)
(primes = restriction to the first n indices), via finSumFinEquiv + det_fromBlocks₂₂ + CauchyDetSchur (S2912).
-/
import Mathlib
import CauchyDetSchur

open Matrix

namespace CauchyDet

variable {K : Type*} [Field K]

def cauchy {n : ℕ} (x y : Fin n → K) : Matrix (Fin n) (Fin n) K := Matrix.of fun i j => (x i + y j)⁻¹

theorem det_cauchy_succ {n : ℕ} (x y : Fin (n + 1) → K) (h : ∀ i j, x i + y j ≠ 0) :
    (cauchy x y).det =
      (x (Fin.last n) + y (Fin.last n))⁻¹ *
        ((∏ i : Fin n, (x (Fin.last n) - x i.castSucc) / (x i.castSucc + y (Fin.last n))) *
          (cauchy (x ∘ Fin.castSucc) (y ∘ Fin.castSucc)).det *
          ∏ j : Fin n, (y (Fin.last n) - y j.castSucc) / (x (Fin.last n) + y j.castSucc)) := by
  set xn := x (Fin.last n); set yn := y (Fin.last n)
  have hd : xn + yn ≠ 0 := h _ _
  -- block form
  let e : Fin n ⊕ Fin 1 ≃ Fin (n + 1) := finSumFinEquiv
  have hblk : (cauchy x y).submatrix e e =
      Matrix.fromBlocks (cauchy (x ∘ Fin.castSucc) (y ∘ Fin.castSucc))
        (Matrix.of fun (i : Fin n) (_ : Fin 1) => (x i.castSucc + yn)⁻¹)
        (Matrix.of fun (_ : Fin 1) (j : Fin n) => (xn + y j.castSucc)⁻¹)
        (Matrix.of fun (_ : Fin 1) (_ : Fin 1) => (xn + yn)⁻¹) := by
    ext (i | i) (j | j) <;> simp [cauchy, e, Fin.fin_one_eq_zero, xn, yn] <;> rfl
  set D : Matrix (Fin 1) (Fin 1) K := Matrix.of fun _ _ => (xn + yn)⁻¹
  have hDdet : D.det = (xn + yn)⁻¹ := by simp [D, Matrix.det_unique]
  letI : Invertible D := Matrix.invertibleOfIsUnitDet D (by rw [hDdet]; exact isUnit_iff_ne_zero.mpr (inv_ne_zero hd))
  have hinv : ⅟D = Matrix.of fun _ _ => (xn + yn) := by
    rw [Matrix.invOf_eq_nonsing_inv]
    apply Matrix.inv_eq_left_inv
    ext i j; fin_cases i; fin_cases j; simp [D, Matrix.mul_apply, hd]
  rw [← Matrix.det_submatrix_equiv_self e, hblk, Matrix.det_fromBlocks₂₂, hDdet, hinv]
  congr 1
  have hschur := schur_matrix (x ∘ Fin.castSucc) (y ∘ Fin.castSucc) xn yn (fun i j => h _ _)
    (fun i => h _ _) (fun j => h _ _) hd
  have hsub : cauchy (x ∘ Fin.castSucc) (y ∘ Fin.castSucc) -
      Matrix.of (fun (i : Fin n) (_ : Fin 1) => (x i.castSucc + yn)⁻¹) * Matrix.of (fun (_ : Fin 1) (_ : Fin 1) => xn + yn) *
        Matrix.of (fun (_ : Fin 1) (j : Fin n) => (xn + y j.castSucc)⁻¹)
      = Matrix.of fun i j => (x i.castSucc + y j.castSucc)⁻¹ - (x i.castSucc + yn)⁻¹ * (xn + yn) * (xn + y j.castSucc)⁻¹ := by
    ext i j; simp [cauchy, Matrix.mul_apply]
  rw [hsub]
  simp only [Function.comp] at hschur
  rw [hschur, Matrix.det_mul, Matrix.det_mul, Matrix.det_diagonal, Matrix.det_diagonal]
  simp [cauchy, Function.comp]

end CauchyDet

#print axioms CauchyDet.det_cauchy_succ
