import HodgeIndexInequalityCore
import HodgeIndexOneOnModel

set_option maxHeartbeats 40000000

open IsDedekindDomain NumberField QuadraticMap Submodule

noncomputable section

variable (K : Type*) [Field K] [NumberField K]

namespace ArakelovVecR

/-!
# §1721 ★★★★★ — THE ARITHMETIC HODGE INDEX INEQUALITY on the Arakelov ℝ-model
On any finite-dim subspace `W` carrying a hodge-positive vector, for any `u, v ∈ W` with the Hodge
self-intersection `hodgeBilin u u > 0`:  `hodgeBilin(u,v)² ≥ hodgeBilin(u,u) · hodgeBilin(v,v)`.
This is the Lorentzian positivity of the arithmetic intersection form — the exact analog of the
inequality that yields RH in the function-field case. Assembled from §1719 (sigPos=1) + §1721a.
-/

theorem hodge_index_inequality_model
    (W : Submodule ℝ (ArakelovVecR K)) [FiniteDimensional ℝ W]
    (wp : ArakelovVecR K) (hwpW : wp ∈ W) (hwppos : 0 < hodgeBilin K wp wp)
    (u v : ArakelovVecR K) (huW : u ∈ W) (hvW : v ∈ W)
    (hupos : 0 < hodgeBilin K u u) :
    hodgeBilin K u u * hodgeBilin K v v ≤ hodgeBilin K u v ^ 2 := by
  set Q := ((hodgeBilin K).restrict W).toQuadraticMap with hQdef
  have hres : ∀ z : W, Q z = hodgeBilin K (z : ArakelovVecR K) (z : ArakelovVecR K) := by
    intro z; rw [hQdef, LinearMap.BilinMap.toQuadraticMap_apply]; rfl
  have hsig : sigPos (R := ℝ) (M := W) Q = 1 := sigPos_restrict_eq_one K W wp hwpW hwppos
  -- polar relates to hodgeBilin (symmetric): polar Q u' v' = 2 * hodgeBilin u v
  set u' : W := ⟨u, huW⟩ with hu'
  set v' : W := ⟨v, hvW⟩ with hv'
  have hpolar : polar Q u' v' = 2 * hodgeBilin K u v := by
    rw [hQdef, LinearMap.BilinMap.polar_toQuadraticMap]
    show ((hodgeBilin K).restrict W) u' v' + ((hodgeBilin K).restrict W) v' u'
      = 2 * hodgeBilin K u v
    have h1 : ((hodgeBilin K).restrict W) u' v' = hodgeBilin K u v := rfl
    have h2 : ((hodgeBilin K).restrict W) v' u' = hodgeBilin K v u := rfl
    rw [h1, h2, hodgeBilin_symm K v u]; ring
  have hQu : Q u' = hodgeBilin K u u := hres u'
  have hQv : Q v' = hodgeBilin K v v := hres v'
  have hQupos : 0 < Q u' := by rw [hQu]; exact hupos
  have hineq := hodge_index_inequality Q hsig u' v' hQupos
  rw [hpolar, hQu, hQv] at hineq
  nlinarith [hineq]

#print axioms hodge_index_inequality_model
