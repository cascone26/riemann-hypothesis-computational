import PSDNullRadical

/-!
# §1747 — For a PSD form, the radical is exactly the kernel of the form-as-map

A positive-semidefinite symmetric bilinear form `B : M →ₗ M →ₗ ℝ` has

`B u = 0  ⟺  B u u = 0`,

i.e. the null cone `{u | B u u = 0}` coincides with the kernel of `B` viewed as the linear map `M → M*`.  (`⟸`
is §1744: a null vector pairs to zero with everything, so its associated functional is zero.)  So `rad(B) = ker B`
as submodules — the degenerate directions of the quadratic form are precisely the vectors the form cannot see.
For the PSD-saturated Weil form this identifies its radical with the kernel of the pairing operator.  EXIT0.
-/

/-- **Radical = kernel.**  For PSD symmetric `B`, the functional `B u` is zero iff `u` is null. -/
theorem psd_ker_eq_null {M : Type*} [AddCommGroup M] [Module ℝ M]
    (B : LinearMap.BilinForm ℝ M) (hsymm : ∀ x y, B x y = B y x)
    (hpsd : ∀ x, 0 ≤ B x x) (u : M) :
    B u = 0 ↔ B u u = 0 := by
  constructor
  · intro h
    rw [h]; rfl
  · intro h
    ext w
    exact psd_null_in_radical B hsymm hpsd u h w
