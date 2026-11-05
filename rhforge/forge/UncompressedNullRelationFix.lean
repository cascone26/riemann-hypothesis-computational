import Eq7p17DetIdentity
import EOffInterpVector
import STOffFinalAssembly
import DLNKOff0

open scoped MatrixOrder Matrix
open Matrix

noncomputable def A_off0_full (K : ℕ) (L : ℝ) (z : Fin K → ℂ) :
    Matrix (Fin (2 * (2 * K) + 1)) (Fin (2 * (2 * K) + 1)) ℝ :=
  fun m n => (2 * Real.pi * (((n : ℕ) : ℝ) - (2 * K : ℝ)) / L) *
             ((if m = n then (1 : ℝ) else 0) - (c_off K L z m))

lemma S_mulVec_c_off_eq_zero (K : ℕ) (L : ℝ) (hL : 0 < L) (z : Fin K → ℂ)
    (hz : ∀ k, (z k).im ≠ 0) (i : Fin (2 * (2 * K) + 1)) :
    ∑ m, (S_T_off_real (2 * K) L K z) i m * (c_off K L z m) = 0 := by
  have hLne : L ≠ 0 := ne_of_gt hL
  have he := S_T_off_real_mulVec_e_off_eq_zero K L hLne z hz i
  rw [Matrix.mulVec, dotProduct] at he
  simp only [e_off] at he
  have hsqrtL : Real.sqrt L ≠ 0 := Real.sqrt_ne_zero'.mpr hL
  have hrw : ∑ m, (S_T_off_real (2 * K) L K z) i m * (Real.sqrt L * c_off K L z m)
      = Real.sqrt L * ∑ m, (S_T_off_real (2 * K) L K z) i m * (c_off K L z m) := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl (fun m _ => by ring)
  rw [hrw] at he
  rcases mul_eq_zero.mp he with h | h
  · exact absurd h hsqrtL
  · exact h

theorem uncompressed_null_relation_core (K : ℕ) (L : ℝ) (hL : 0 < L) (z : Fin K → ℂ)
    (hz : ∀ k, (z k).im ≠ 0) :
    S_T_off_real (2 * K) L K z * (A_off0_full K L z) =
    S_T_off_real (2 * K) L K z * (D_LN (2 * K) L) := by
  have hcz := S_mulVec_c_off_eq_zero K L hL z hz
  ext i j
  rw [Matrix.mul_apply, Matrix.mul_apply]
  have hD : ∀ m : Fin (2 * (2 * K) + 1),
      (D_LN (2 * K) L) m j = if m = j then (D_LN (2 * K) L) j j else 0 := by
    intro m
    by_cases h : m = j
    · simp [h]
    · unfold D_LN
      simp [Matrix.diagonal_apply, h]
  have hRHS : ∑ m, (S_T_off_real (2 * K) L K z) i m * (D_LN (2 * K) L) m j
      = (S_T_off_real (2 * K) L K z) i j * (D_LN (2 * K) L) j j := by
    calc ∑ m, (S_T_off_real (2 * K) L K z) i m * (D_LN (2 * K) L) m j
        = ∑ m, (S_T_off_real (2 * K) L K z) i m *
            (if m = j then (D_LN (2 * K) L) j j else 0) :=
          Finset.sum_congr rfl (fun m _ => by rw [hD m])
      _ = ∑ m, (if m = j then (S_T_off_real (2 * K) L K z) i m * (D_LN (2 * K) L) j j else 0) :=
          Finset.sum_congr rfl (fun m _ => by split_ifs <;> ring)
      _ = (S_T_off_real (2 * K) L K z) i j * (D_LN (2 * K) L) j j := by
          simp [Finset.sum_ite_eq']
  rw [hRHS]
  set c : ℝ := 2 * Real.pi * (((j : ℕ) : ℝ) - (2 * K : ℝ)) / L with hc_def
  have hDjj : (D_LN (2 * K) L) j j = c := by
    unfold D_LN
    simp [Matrix.diagonal_apply, hc_def]
  rw [hDjj]
  have key : ∀ m : Fin (2 * (2 * K) + 1),
      (S_T_off_real (2 * K) L K z) i m * (A_off0_full K L z m j)
      = (if m = j then c * (S_T_off_real (2 * K) L K z) i m else 0)
        - c * ((S_T_off_real (2 * K) L K z) i m * c_off K L z m) := by
    intro m
    show (S_T_off_real (2 * K) L K z) i m *
      ((2 * Real.pi * (((j : ℕ) : ℝ) - (2 * K : ℝ)) / L) *
        ((if m = j then (1 : ℝ) else 0) - c_off K L z m)) = _
    rw [← hc_def]
    split_ifs <;> ring
  rw [Finset.sum_congr rfl (fun m _ => key m)]
  rw [Finset.sum_sub_distrib]
  rw [Finset.sum_ite_eq' Finset.univ j (fun m => c * (S_T_off_real (2 * K) L K z) i m)]
  simp only [Finset.mem_univ, if_true]
  rw [← Finset.mul_sum, hcz i]
  ring

#print axioms uncompressed_null_relation_core
