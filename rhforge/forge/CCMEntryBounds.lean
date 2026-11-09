import CCMWeightGram

/-!
# S2948 — entrywise structure of the CCM zero-side Gram (lane 2)

For G = C·diag(w)·Cᵀ with w ≥ 0 (in particular the CCM/Shi weights):
* `gram_entry` : G_mn = Σ_j w_j C_mj C_nj;
* `gram_diag_nonneg` : G_mm ≥ 0;
* `gram_entry_sq_le` : G_mn² ≤ G_mm·G_nn  (weighted Cauchy–Schwarz — the 2×2-minor form of PSD);
* `ccm_diag_le` : G^ccm_mm ≤ (4/L)·Σ_j C_mj².
For the Shi/CCM Cauchy matrix C_mj = 1/(ν_m+γ_j) these say the zero-side coupling between test functions m, n is
controlled by their self-energies.
-/

open Matrix Finset

namespace CCMEntry

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq n] [DecidableEq m]

theorem gram_entry (C : Matrix n m ℝ) (w : m → ℝ) (a b : n) :
    (C * diagonal w * Cᵀ) a b = ∑ j, w j * C a j * C b j := by
  simp only [mul_apply, transpose_apply]
  refine sum_congr rfl fun j _ => ?_
  simp [diagonal_apply]; try ring_nf; try simp

theorem gram_diag_nonneg (C : Matrix n m ℝ) {w : m → ℝ} (hw : ∀ j, 0 ≤ w j) (a : n) :
    0 ≤ (C * diagonal w * Cᵀ) a a := by
  rw [gram_entry]
  exact sum_nonneg fun j _ => by nlinarith [hw j, mul_self_nonneg (C a j)]

theorem gram_entry_sq_le (C : Matrix n m ℝ) {w : m → ℝ} (hw : ∀ j, 0 ≤ w j) (a b : n) :
    ((C * diagonal w * Cᵀ) a b) ^ 2
      ≤ (C * diagonal w * Cᵀ) a a * (C * diagonal w * Cᵀ) b b := by
  rw [gram_entry, gram_entry, gram_entry]
  have h := sum_mul_sq_le_sq_mul_sq univ (fun j => Real.sqrt (w j) * C a j)
    (fun j => Real.sqrt (w j) * C b j)
  have hs : ∀ j, Real.sqrt (w j) ^ 2 = w j := fun j => Real.sq_sqrt (hw j)
  have e1 : ∑ j, Real.sqrt (w j) * C a j * (Real.sqrt (w j) * C b j) = ∑ j, w j * C a j * C b j :=
    sum_congr rfl fun j _ => by linear_combination (C a j * C b j) * hs j
  have e2 : ∀ x : n, ∑ j, (Real.sqrt (w j) * C x j) ^ 2 = ∑ j, w j * C x j * C x j :=
    fun x => sum_congr rfl fun j _ => by linear_combination (C x j ^ 2) * hs j
  rw [e1, e2, e2] at h
  exact h

theorem ccm_diag_le {L : ℝ} (hL : 0 < L) (C : Matrix n m ℝ) (γ : m → ℝ) (a : n) :
    (C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ) a a ≤ 4 / L * ∑ j, C a j ^ 2 := by
  rw [gram_entry, mul_sum]
  exact sum_le_sum fun j _ => by
    have := CCMWeight.ccmW_le hL (γ j)
    nlinarith [sq_nonneg (C a j)]

theorem ccm_entry_sq_le {L : ℝ} (hL : 0 < L) (C : Matrix n m ℝ) (γ : m → ℝ) (a b : n) :
    ((C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ) a b) ^ 2
      ≤ (C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ) a a
        * (C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ) b b :=
  gram_entry_sq_le C (fun j => CCMWeight.ccmW_nonneg hL (γ j)) a b

end CCMEntry

#print axioms CCMEntry.gram_entry_sq_le
#print axioms CCMEntry.ccm_diag_le
#print axioms CCMEntry.ccm_entry_sq_le
