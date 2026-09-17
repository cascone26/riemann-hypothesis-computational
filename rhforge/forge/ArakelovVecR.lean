import Mathlib.RingTheory.DedekindDomain.AdicValuation
import Mathlib.NumberTheory.NumberField.InfinitePlace.Basic
import Mathlib.LinearAlgebra.BilinearMap
import Mathlib.Tactic

set_option maxHeartbeats 1000000

open IsDedekindDomain NumberField

noncomputable section

variable (K : Type*) [Field K] [NumberField K]

/-!
# §1712 — THE ℝ-BASE-CHANGED ARAKELOV SPACE: a genuine `Module ℝ` with the pairing as a `BilinForm`

The §1659 Hodge-signature theorem consumes a `LinearMap.BilinForm ℝ M` on a real module; the
concrete divisor model (§1652/§1710) is a ℤ-Finsupp × ℝ-functions hybrid with no scalar action.
This brick base-changes: `ArakelovVecR := (HeightOneSpectrum 𝒪_K →₀ ℝ) × (InfinitePlace K → ℝ)`
(module instances for free), and packages the intersection pairing as an honest bilinear form
`pairingBilin : LinearMap.BilinForm ℝ (ArakelovVecR K)`, symmetric.  The finite-dimensional
subspaces of this space are where §1659 will be instantiated (scoped as the next focused build).
-/

abbrev ArakelovVecR :=
  (HeightOneSpectrum (RingOfIntegers K) →₀ ℝ) × (InfinitePlace K → ℝ)

namespace ArakelovVecR

/-- The pairing as a raw function: ℓ²-type on both parts, multiplicity-weighted at infinity. -/
noncomputable def pairingFun (D E : ArakelovVecR K) : ℝ :=
  D.1.sum (fun v a => a * E.1 v) + ∑ w : InfinitePlace K, D.2 w * E.2 w * (w.mult : ℝ)

theorem pairingFun_add_left (D D' E : ArakelovVecR K) :
    pairingFun K (D + D') E = pairingFun K D E + pairingFun K D' E := by
  unfold pairingFun
  have hfin : (D.1 + D'.1).sum (fun v a => a * E.1 v)
      = D.1.sum (fun v a => a * E.1 v) + D'.1.sum (fun v a => a * E.1 v) :=
    Finsupp.sum_add_index' (fun v => zero_mul _) (fun v b₁ b₂ => add_mul _ _ _)
  have h1 : (D + D').1 = D.1 + D'.1 := rfl
  have h2 : ∀ w, (D + D').2 w = D.2 w + D'.2 w := fun _ => rfl
  rw [h1, hfin]
  have hinf : (∑ w : InfinitePlace K, (D + D').2 w * E.2 w * (w.mult : ℝ))
      = (∑ w : InfinitePlace K, D.2 w * E.2 w * (w.mult : ℝ))
        + ∑ w : InfinitePlace K, D'.2 w * E.2 w * (w.mult : ℝ) := by
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun w _ => by rw [h2]; ring
  rw [hinf]
  ring

theorem pairingFun_smul_left (c : ℝ) (D E : ArakelovVecR K) :
    pairingFun K (c • D) E = c * pairingFun K D E := by
  unfold pairingFun
  have h1 : (c • D).1 = c • D.1 := rfl
  have h2 : ∀ w, (c • D).2 w = c * D.2 w := fun _ => rfl
  have hfin : (c • D.1).sum (fun v a => a * E.1 v)
      = c * D.1.sum (fun v a => a * E.1 v) := by
    rw [Finsupp.sum_smul_index (fun v => zero_mul _), Finsupp.mul_sum]
    exact Finsupp.sum_congr fun v _ => by ring
  rw [h1, hfin]
  have hinf : (∑ w : InfinitePlace K, (c • D).2 w * E.2 w * (w.mult : ℝ))
      = c * ∑ w : InfinitePlace K, D.2 w * E.2 w * (w.mult : ℝ) := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun w _ => by rw [h2]; ring
  rw [hinf]
  ring

theorem pairingFun_symm (D E : ArakelovVecR K) :
    pairingFun K D E = pairingFun K E D := by
  unfold pairingFun
  classical
  have hfin : D.1.sum (fun v a => a * E.1 v) = E.1.sum (fun v a => a * D.1 v) := by
    rw [Finsupp.sum_of_support_subset D.1
          (Finset.subset_union_left : D.1.support ⊆ D.1.support ∪ E.1.support)
          (fun v a => a * E.1 v) (fun v _ => by ring),
        Finsupp.sum_of_support_subset E.1
          (Finset.subset_union_right : E.1.support ⊆ D.1.support ∪ E.1.support)
          (fun v a => a * D.1 v) (fun v _ => by ring)]
    exact Finset.sum_congr rfl fun v _ => by ring
  rw [hfin]
  have hinf : (∑ w : InfinitePlace K, D.2 w * E.2 w * (w.mult : ℝ))
      = ∑ w : InfinitePlace K, E.2 w * D.2 w * (w.mult : ℝ) :=
    Finset.sum_congr rfl fun w _ => by ring
  rw [hinf]

theorem pairingFun_add_right (D E E' : ArakelovVecR K) :
    pairingFun K D (E + E') = pairingFun K D E + pairingFun K D E' := by
  rw [pairingFun_symm, pairingFun_add_left, pairingFun_symm K E D, pairingFun_symm K E' D]

theorem pairingFun_smul_right (c : ℝ) (D E : ArakelovVecR K) :
    pairingFun K D (c • E) = c * pairingFun K D E := by
  rw [pairingFun_symm, pairingFun_smul_left, pairingFun_symm K E D]

/-- The intersection pairing as an honest bilinear form on the ℝ-model. -/
noncomputable def pairingBilin : LinearMap.BilinForm ℝ (ArakelovVecR K) :=
  LinearMap.mk₂ ℝ (pairingFun K)
    (pairingFun_add_left K)
    (fun c D E => by rw [pairingFun_smul_left]; rfl)
    (pairingFun_add_right K)
    (fun c D E => by rw [pairingFun_smul_right]; rfl)

theorem pairingBilin_apply (D E : ArakelovVecR K) :
    pairingBilin K D E = pairingFun K D E := rfl

theorem pairingBilin_symm (D E : ArakelovVecR K) :
    pairingBilin K D E = pairingBilin K E D := pairingFun_symm K D E

#print axioms pairingBilin_symm
