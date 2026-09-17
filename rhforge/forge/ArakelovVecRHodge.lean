import ArakelovVecR
import Mathlib.RingTheory.Ideal.Norm.AbsNorm

set_option maxHeartbeats 1000000

open IsDedekindDomain NumberField

noncomputable section

variable (K : Type*) [Field K] [NumberField K]

namespace ArakelovVecR

/-!
# §1713 — DEG FUNCTIONAL + POSDEF + THE HODGE FORM on the ℝ-model

Completes the ℝ-model package: `degLin` (the full Arakelov degree as a genuine `→ₗ[ℝ] ℝ`),
positive-definiteness of `pairingBilin` (ℓ²-type), and the Hodge form
`hodgeBilin = deg ⊗ deg − pairing` as a symmetric `BilinForm` that is **strictly negative on the
nonzero kernel of `degLin`** — §1659's negdef-primitive hypothesis on the base-changed space.
Remaining for the focused instantiation session: choose a finite-dimensional subspace containing a
hyperbolic pair and apply `hodge_signature_of_adjunction`.
-/

/-- The full Arakelov degree as a raw function on the ℝ-model. -/
noncomputable def degFun (D : ArakelovVecR K) : ℝ :=
  D.1.sum (fun v a => a * Real.log (Ideal.absNorm v.asIdeal : ℝ))
    + ∑ w : InfinitePlace K, (w.mult : ℝ) * D.2 w

theorem degFun_add (D E : ArakelovVecR K) :
    degFun K (D + E) = degFun K D + degFun K E := by
  unfold degFun
  have hfin : (D.1 + E.1).sum (fun v a => a * Real.log (Ideal.absNorm v.asIdeal : ℝ))
      = D.1.sum (fun v a => a * Real.log (Ideal.absNorm v.asIdeal : ℝ))
        + E.1.sum (fun v a => a * Real.log (Ideal.absNorm v.asIdeal : ℝ)) :=
    Finsupp.sum_add_index' (fun v => zero_mul _) (fun v b₁ b₂ => add_mul _ _ _)
  have h1 : (D + E).1 = D.1 + E.1 := rfl
  have h2 : ∀ w, (D + E).2 w = D.2 w + E.2 w := fun _ => rfl
  rw [h1, hfin]
  have hinf : (∑ w : InfinitePlace K, (w.mult : ℝ) * (D + E).2 w)
      = (∑ w : InfinitePlace K, (w.mult : ℝ) * D.2 w)
        + ∑ w : InfinitePlace K, (w.mult : ℝ) * E.2 w := by
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun w _ => by rw [h2]; ring
  rw [hinf]
  ring

theorem degFun_smul (c : ℝ) (D : ArakelovVecR K) :
    degFun K (c • D) = c * degFun K D := by
  unfold degFun
  have h1 : (c • D).1 = c • D.1 := rfl
  have h2 : ∀ w, (c • D).2 w = c * D.2 w := fun _ => rfl
  have hfin : (c • D.1).sum (fun v a => a * Real.log (Ideal.absNorm v.asIdeal : ℝ))
      = c * D.1.sum (fun v a => a * Real.log (Ideal.absNorm v.asIdeal : ℝ)) := by
    rw [Finsupp.sum_smul_index (fun v => zero_mul _), Finsupp.mul_sum]
    exact Finsupp.sum_congr fun v _ => by ring
  rw [h1, hfin]
  have hinf : (∑ w : InfinitePlace K, (w.mult : ℝ) * (c • D).2 w)
      = c * ∑ w : InfinitePlace K, (w.mult : ℝ) * D.2 w := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun w _ => by rw [h2]; ring
  rw [hinf]
  ring

/-- The full Arakelov degree as a linear functional. -/
noncomputable def degLin : ArakelovVecR K →ₗ[ℝ] ℝ where
  toFun := degFun K
  map_add' := degFun_add K
  map_smul' := fun c D => by rw [degFun_smul]; rfl

theorem degLin_apply (D : ArakelovVecR K) : degLin K D = degFun K D := rfl

/-- Positive semidefiniteness of the pairing. -/
theorem pairingFun_self_nonneg (D : ArakelovVecR K) : 0 ≤ pairingFun K D D := by
  unfold pairingFun
  apply add_nonneg
  · rw [Finsupp.sum]
    exact Finset.sum_nonneg fun v _ => mul_self_nonneg _
  · exact Finset.sum_nonneg fun w _ =>
      mul_nonneg (mul_self_nonneg _) (Nat.cast_nonneg _)

/-- Definiteness: the pairing vanishes on the diagonal only at zero. -/
theorem pairingFun_self_eq_zero_iff (D : ArakelovVecR K) :
    pairingFun K D D = 0 ↔ D = 0 := by
  constructor
  · intro h
    unfold pairingFun at h
    have hfin0 : (0:ℝ) ≤ D.1.sum (fun v a => a * D.1 v) := by
      rw [Finsupp.sum]
      exact Finset.sum_nonneg fun v _ => mul_self_nonneg _
    have hinf0 : (0:ℝ) ≤ ∑ w : InfinitePlace K, D.2 w * D.2 w * (w.mult : ℝ) :=
      Finset.sum_nonneg fun w _ => mul_nonneg (mul_self_nonneg _) (Nat.cast_nonneg _)
    have hf : D.1.sum (fun v a => a * D.1 v) = 0 := by linarith
    have hi : (∑ w : InfinitePlace K, D.2 w * D.2 w * (w.mult : ℝ)) = 0 := by linarith
    have hD1 : D.1 = 0 := by
      rw [Finsupp.sum] at hf
      have hall := (Finset.sum_eq_zero_iff_of_nonneg
        (fun v _ => mul_self_nonneg (D.1 v))).mp hf
      apply Finsupp.support_eq_empty.mp
      apply Finset.eq_empty_of_forall_notMem
      intro v hv
      have h2 : D.1 v = 0 := mul_self_eq_zero.mp (hall v hv)
      exact (Finsupp.mem_support_iff.mp hv) h2
    have hD2 : D.2 = 0 := by
      have hall := (Finset.sum_eq_zero_iff_of_nonneg
        (fun w _ => mul_nonneg (mul_self_nonneg (D.2 w)) (Nat.cast_nonneg w.mult))).mp hi
      funext w
      have h1 := hall w (Finset.mem_univ w)
      have hm : (0:ℝ) < (w.mult : ℝ) := by exact_mod_cast InfinitePlace.mult_pos
      have h2 : D.2 w * D.2 w = 0 := by
        by_contra hne
        exact (mul_ne_zero hne hm.ne') h1
      exact mul_self_eq_zero.mp h2
    exact Prod.ext hD1 hD2
  · intro h
    rw [h]
    unfold pairingFun
    simp

theorem pairingFun_self_pos_of_ne_zero {D : ArakelovVecR K} (hD : D ≠ 0) :
    0 < pairingFun K D D :=
  lt_of_le_of_ne (pairingFun_self_nonneg K D)
    (fun h => hD ((pairingFun_self_eq_zero_iff K D).mp h.symm))

/-- The Hodge bilinear form `deg ⊗ deg − pairing` on the ℝ-model. -/
noncomputable def hodgeBilin : LinearMap.BilinForm ℝ (ArakelovVecR K) :=
  LinearMap.mk₂ ℝ (fun D E => degFun K D * degFun K E - pairingFun K D E)
    (fun D D' E => by rw [degFun_add, pairingFun_add_left]; ring)
    (fun c D E => by rw [degFun_smul, pairingFun_smul_left]; simp; ring)
    (fun D E E' => by rw [degFun_add, pairingFun_add_right]; ring)
    (fun c D E => by rw [degFun_smul, pairingFun_smul_right]; simp; ring)

theorem hodgeBilin_apply (D E : ArakelovVecR K) :
    hodgeBilin K D E = degFun K D * degFun K E - pairingFun K D E := rfl

theorem hodgeBilin_symm (D E : ArakelovVecR K) :
    hodgeBilin K D E = hodgeBilin K E D := by
  rw [hodgeBilin_apply, hodgeBilin_apply, pairingFun_symm]
  ring

/-- **Negdef on the primitive part**: the Hodge form is strictly negative on every nonzero
element of `ker degLin` — §1659's hypothesis on the base-changed model. -/
theorem hodgeBilin_self_neg_of_deg_zero {D : ArakelovVecR K}
    (hdeg : degLin K D = 0) (hD : D ≠ 0) :
    hodgeBilin K D D < 0 := by
  rw [hodgeBilin_apply]
  rw [degLin_apply] at hdeg
  rw [hdeg]
  have := pairingFun_self_pos_of_ne_zero K hD
  linarith

#print axioms hodgeBilin_self_neg_of_deg_zero
#print axioms pairingFun_self_eq_zero_iff
