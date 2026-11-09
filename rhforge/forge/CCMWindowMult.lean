import CCMResonanceGap

/-!
# S2964 — window multiplication: G_ccm(kL) ⪯ k·G_ccm(L) (lane 2)

`|sin(kx)| ≤ k|sin x|` (induction on k via the addition formula) gives, for the CCM zero-side weight
`w_L(γ) = 4 sin²(γL/2)/L`, the pointwise inequality `w_{kL}(γ) ≤ k·w_L(γ)`.  Hence:
* `abs_sin_nat_mul_le` — `|sin(k·x)| ≤ k·|sin x|`.
* `ccmW_nat_mul_le` — `ccmW (k·L) γ ≤ k·ccmW L γ`  (k ≥ 1, L > 0).
* `ccm_gram_nat_mul_le` — `k·G_ccm(L) − G_ccm(kL) ⪰ 0`.
* `ccm_rayleigh_nat_mul_le` — Rayleigh form: xᵀG(kL)x ≤ k·xᵀG(L)x.
Reading: zero-side margins can shrink under window multiplication only by a factor ≤ k, i.e. any margin decay
faster than linear in the window must come from the arithmetic side (W02 − W_R − W_p), not from the zero side.
-/

open Matrix Finset Real

namespace CCMWindow

theorem abs_sin_nat_mul_le (k : ℕ) (x : ℝ) : |Real.sin (k * x)| ≤ k * |Real.sin x| := by
  induction k with
  | zero => simp
  | succ k ih =>
    have hk : ((k + 1 : ℕ) : ℝ) * x = (k : ℝ) * x + x := by push_cast; ring
    rw [hk, Real.sin_add]
    have h1 := Real.abs_cos_le_one x
    have h2 := Real.abs_cos_le_one ((k : ℝ) * x)
    have e1 := mul_le_mul_of_nonneg_left h1 (abs_nonneg (Real.sin ((k : ℝ) * x)))
    have e2 := mul_le_mul_of_nonneg_right h2 (abs_nonneg (Real.sin x))
    calc |Real.sin ((k : ℝ) * x) * Real.cos x + Real.cos ((k : ℝ) * x) * Real.sin x|
        ≤ |Real.sin ((k : ℝ) * x) * Real.cos x| + |Real.cos ((k : ℝ) * x) * Real.sin x| :=
          abs_add_le _ _
      _ ≤ |Real.sin ((k : ℝ) * x)| + |Real.sin x| := by
          rw [abs_mul, abs_mul]; linarith
      _ ≤ ((k + 1 : ℕ) : ℝ) * |Real.sin x| := by push_cast; linarith

theorem ccmW_nat_mul_le {L : ℝ} (hL : 0 < L) {k : ℕ} (hk : 1 ≤ k) (γ : ℝ) :
    CCMWeight.ccmW ((k : ℝ) * L) γ ≤ (k : ℝ) * CCMWeight.ccmW L γ := by
  unfold CCMWeight.ccmW
  have hkR : (1 : ℝ) ≤ (k : ℝ) := by exact_mod_cast hk
  have hkpos : (0 : ℝ) < (k : ℝ) := by linarith
  set y : ℝ := γ * L / 2 with hy
  rw [show γ * ((k : ℝ) * L) / 2 = (k : ℝ) * y by rw [hy]; ring]
  have h := abs_sin_nat_mul_le k y
  have hsq : Real.sin ((k : ℝ) * y) ^ 2 ≤ (k : ℝ) ^ 2 * Real.sin y ^ 2 := by
    have := pow_le_pow_left₀ (abs_nonneg _) h 2
    rw [sq_abs, mul_pow, sq_abs] at this
    exact this
  have hkL : 0 < (k : ℝ) * L := by positivity
  calc 4 * Real.sin ((k : ℝ) * y) ^ 2 / ((k : ℝ) * L)
      ≤ 4 * ((k : ℝ) ^ 2 * Real.sin y ^ 2) / ((k : ℝ) * L) := by gcongr
    _ = (k : ℝ) * (4 * Real.sin y ^ 2 / L) := by field_simp

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq n] [DecidableEq m]

theorem ccm_gram_nat_mul_le {L : ℝ} (hL : 0 < L) {k : ℕ} (hk : 1 ≤ k) (C : Matrix n m ℝ)
    (γ : m → ℝ) :
    ((k : ℝ) • (C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ)
      - C * diagonal (fun j => CCMWeight.ccmW ((k : ℝ) * L) (γ j)) * Cᵀ).PosSemidef := by
  rw [← CCMFejer.gram_smul]
  exact CCMLoewner.gram_mono C (fun j => by
    simpa [Pi.smul_apply, smul_eq_mul] using ccmW_nat_mul_le hL hk (γ j))

theorem ccm_rayleigh_nat_mul_le {L : ℝ} (hL : 0 < L) {k : ℕ} (hk : 1 ≤ k) (C : Matrix n m ℝ)
    (γ : m → ℝ) (x : n → ℝ) :
    x ⬝ᵥ ((C * diagonal (fun j => CCMWeight.ccmW ((k : ℝ) * L) (γ j)) * Cᵀ) *ᵥ x)
      ≤ (k : ℝ) * (x ⬝ᵥ ((C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ) *ᵥ x)) := by
  have h := CCMResonance.rayleigh_mono C (w' := (k : ℝ) • fun j => CCMWeight.ccmW L (γ j))
    (fun j => by simpa [Pi.smul_apply, smul_eq_mul] using ccmW_nat_mul_le hL hk (γ j)) x
  rw [CCMFejer.gram_smul, Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul] at h
  exact h

end CCMWindow

#print axioms CCMWindow.abs_sin_nat_mul_le
#print axioms CCMWindow.ccmW_nat_mul_le
#print axioms CCMWindow.ccm_gram_nat_mul_le
#print axioms CCMWindow.ccm_rayleigh_nat_mul_le
