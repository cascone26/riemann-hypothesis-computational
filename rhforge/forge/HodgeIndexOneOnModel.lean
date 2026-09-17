import UniformHodgeIndexCore
import ArakelovVecRHodge

set_option maxHeartbeats 40000000

open IsDedekindDomain NumberField QuadraticMap Submodule

noncomputable section

variable (K : Type*) [Field K] [NumberField K]

namespace ArakelovVecR

/-!
# §1719 ★★★★ — ARITHMETIC HODGE INDEX ON EVERY SUBSPACE (ℝ-model): `sigPos = 1`
Delegates the sigPos linear algebra to abstract §1719a (general ordered field, plain-inequality
negativity hypothesis), so no ℝ-instance diamond / Neg-instance issue arises on this side.
-/

theorem sigPos_restrict_eq_one
    (W : Submodule ℝ (ArakelovVecR K)) [FiniteDimensional ℝ W]
    (w : ArakelovVecR K) (hwW : w ∈ W) (hwpos : 0 < hodgeBilin K w w) :
    sigPos (R := ℝ) (M := W) (((hodgeBilin K).restrict W).toQuadraticMap) = 1 := by
  set Q := ((hodgeBilin K).restrict W).toQuadraticMap with hQdef
  set φ : W →ₗ[ℝ] ℝ := (degLin K).comp W.subtype with hφdef
  have hres : ∀ z : W, Q z = hodgeBilin K (z : ArakelovVecR K) (z : ArakelovVecR K) := by
    intro z
    rw [hQdef, LinearMap.BilinMap.toQuadraticMap_apply]
    rfl
  have hφap : ∀ z : W, φ z = degLin K (z : ArakelovVecR K) := fun z => rfl
  have hneg : ∀ z ∈ LinearMap.ker φ, z ≠ 0 → Q z < 0 := by
    intro z hzmem hz0
    have hz0' : ((z : W) : ArakelovVecR K) ≠ 0 := fun h => hz0 (Subtype.ext h)
    have hdeg0 : degLin K (((z : W) : ArakelovVecR K)) = 0 := by
      have := LinearMap.mem_ker.mp hzmem; rwa [hφap] at this
    rw [hres]
    exact hodgeBilin_self_neg_of_deg_zero K hdeg0 hz0'
  have hQwv : 0 < Q (⟨w, hwW⟩ : W) := by rw [hres]; exact hwpos
  have hwv0 : (⟨w, hwW⟩ : W) ≠ 0 := by
    intro h; rw [h, QuadraticMap.map_zero] at hQwv; exact lt_irrefl 0 hQwv
  have hφwv : φ (⟨w, hwW⟩ : W) ≠ 0 := by
    intro h0
    have hmem : (⟨w, hwW⟩ : W) ∈ LinearMap.ker φ := LinearMap.mem_ker.mpr h0
    have := hneg _ hmem hwv0
    linarith [hQwv]
  exact sigPos_eq_one_of_pos_and_negdef_ker Q φ (⟨w, hwW⟩ : W) hQwv hφwv hneg

#print axioms sigPos_restrict_eq_one
