import PSDKerEqNull

/-!
# §1750 — The descent to `M / rad`: well-definedness + anisotropy

The two ingredients that make "the PSD form descends to a positive-definite form on `M / rad(B)`" precise,
stated at the element level (no quotient-object bookkeeping needed):

* **well-definedness** — `B` is constant along radical cosets in each argument: `u − u' ∈ rad ⟹ B u v = B u' v`;
* **anisotropy** — off the radical the form is strictly positive: `u ∉ rad ⟹ 0 < B u u`.

Together these say the induced form on the free quotient `M / rad(B)` is a genuine positive-definite inner
product — the abstract statement `NeronTateRadical` proves concretely for the Néron–Tate height (anisotropic on
`E(K)/E(K)_tors`).  For the Weil form: the `Λ=0` saturation is entirely inside `rad`, and modulo it the
positivity is strict.  EXIT0 / std-axioms.
-/

variable {M : Type*} [AddCommGroup M] [Module ℝ M]
  (B : LinearMap.BilinForm ℝ M) (hsymm : ∀ x y, B x y = B y x) (hpsd : ∀ x, 0 ≤ B x x)
  (N : Submodule ℝ M) (hrad : ∀ u, u ∈ N ↔ B u u = 0)

include hsymm hpsd hrad

/-- **Descent well-definedness.**  `B` is constant along radical cosets in the first argument. -/
theorem psd_descent_welldef (u u' v : M) (h : u - u' ∈ N) : B u v = B u' v := by
  have hz : B (u - u') = 0 := (psd_ker_eq_null B hsymm hpsd (u - u')).mpr ((hrad (u - u')).mp h)
  have h0 : B (u - u') v = 0 := by rw [hz]; rfl
  have hlin : B (u - u') v = B u v - B u' v := by simp only [map_sub, LinearMap.sub_apply]
  linarith [h0, hlin]

/-- **Descent anisotropy.**  Off the radical the PSD form is strictly positive. -/
theorem psd_descent_anisotropic (u : M) (hu : u ∉ N) : 0 < B u u :=
  lt_of_le_of_ne (hpsd u) (Ne.symm (fun h => hu ((hrad u).mpr h)))
