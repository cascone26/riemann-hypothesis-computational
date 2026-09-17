import PSDNullRadical

/-!
# §1746 — A PSD form is constant along radical cosets

The descent heart: adding a radical (null) vector never changes the self-energy.

`B` symmetric PSD, `B u u = 0`  ⟹  `∀ w, B (u + w) (u + w) = B w w`.

So `B` is well-defined on cosets `w + rad(B)` — it descends to the quotient `M / rad(B)`, where (by §1743/§1744)
it is positive definite.  For the PSD-saturated Weil form this is the exact statement that only the transverse
(mod-radical) part carries energy — the `Λ=0` saturation is entirely along the radical and invisible to the
descended positive-definite form.  EXIT0 / std-axioms.
-/

/-- **Radical-coset invariance.**  Translating by a null vector preserves the quadratic form. -/
theorem psd_radical_translation_invariant {M : Type*} [AddCommGroup M] [Module ℝ M]
    (B : LinearMap.BilinForm ℝ M) (hsymm : ∀ x y, B x y = B y x)
    (hpsd : ∀ x, 0 ≤ B x x) (u w : M) (hu : B u u = 0) :
    B (u + w) (u + w) = B w w := by
  have h1 : B u w = 0 := psd_null_in_radical B hsymm hpsd u hu w
  have h2 : B w u = 0 := by rw [hsymm w u]; exact h1
  simp only [map_add, LinearMap.add_apply, hu, h1, h2, zero_add, add_zero]
