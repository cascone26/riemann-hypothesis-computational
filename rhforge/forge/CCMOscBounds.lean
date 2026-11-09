import CCMFejerSplit

/-!
# S2953 — Loewner bounds on the OSCILLATORY part of the CCM zero side (lane 2)

With G_cos := C·diag(cos γ_jL)·Cᵀ (S2951: G_ccm = (2/L)(CCᵀ − G_cos)):
* `osc_le` : G_cos ⪯ C·Cᵀ          (cos ≤ 1);
* `osc_ge` : −C·Cᵀ ⪯ G_cos          (cos ≥ −1);
* `osc_rayleigh_abs` : |xᵀG_cos x| ≤ xᵀ(CCᵀ)x.
So the oscillatory Gram is a contraction relative to the fixed Cauchy Gram; positivity of the zero side is
exactly G_cos ⪯ CCᵀ, and the margin is (2/L)·(CCᵀ − G_cos).
-/

open Matrix Finset

namespace CCMOsc

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq n] [DecidableEq m]

theorem osc_le (C : Matrix n m ℝ) (γ : m → ℝ) (L : ℝ) :
    (C * Cᵀ - C * diagonal (fun j => Real.cos (γ j * L)) * Cᵀ).PosSemidef := by
  have h := CCMLoewner.gram_mono C (w := fun j => Real.cos (γ j * L)) (w' := fun _ => (1 : ℝ))
    (fun j => Real.cos_le_one _)
  rwa [CCMLoewner.const_gram, one_smul] at h

theorem osc_ge (C : Matrix n m ℝ) (γ : m → ℝ) (L : ℝ) :
    (C * diagonal (fun j => Real.cos (γ j * L)) * Cᵀ + C * Cᵀ).PosSemidef := by
  have h := CCMLoewner.gram_mono C (w := fun _ => (-1 : ℝ)) (w' := fun j => Real.cos (γ j * L))
    (fun j => Real.neg_one_le_cos _)
  rwa [CCMLoewner.const_gram, neg_one_smul, sub_neg_eq_add] at h

theorem osc_rayleigh_abs (C : Matrix n m ℝ) (γ : m → ℝ) (L : ℝ) (x : n → ℝ) :
    |x ⬝ᵥ ((C * diagonal (fun j => Real.cos (γ j * L)) * Cᵀ) *ᵥ x)| ≤ x ⬝ᵥ ((C * Cᵀ) *ᵥ x) := by
  have h1 := (osc_le C γ L).dotProduct_mulVec_nonneg x
  have h2 := (osc_ge C γ L).dotProduct_mulVec_nonneg x
  simp only [star_trivial, Matrix.sub_mulVec, Matrix.add_mulVec, dotProduct_sub, dotProduct_add] at h1 h2
  rw [abs_le]; constructor <;> linarith

end CCMOsc

#print axioms CCMOsc.osc_le
#print axioms CCMOsc.osc_ge
#print axioms CCMOsc.osc_rayleigh_abs
