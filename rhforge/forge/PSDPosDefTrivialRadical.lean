import PSDNullRadical

/-!
# §1748 — Strict positivity ⟺ trivial radical (the abstract Λ=0 knife-edge)

The clean characterization the saturation phenomenon reduces to: a positive-semidefinite symmetric form is
positive *definite* exactly when its radical is trivial.

`B` PSD symmetric, `(∀ u, B u u = 0 → u = 0)`  ⟹  `∀ u ≠ 0, 0 < B u u`.

So for the Weil form `arch − prime`: strict Weil positivity fails precisely when the radical is nontrivial —
which is the `Λ=0` saturation (verified: PSD with a growing null space).  The de Bruijn–Newman knife-edge
`Λ=0` is, abstractly, exactly the boundary between "radical trivial ⇒ positive definite" and "radical
nontrivial ⇒ merely semidefinite".  EXIT0 / std-axioms.
-/

/-- **Definite from trivial radical.**  A PSD symmetric form with no nonzero null vectors is positive definite. -/
theorem psd_pos_def_of_trivial_radical {M : Type*} [AddCommGroup M] [Module ℝ M]
    (B : LinearMap.BilinForm ℝ M) (hpsd : ∀ x, 0 ≤ B x x)
    (htrivial : ∀ u, B u u = 0 → u = 0) (u : M) (hu : u ≠ 0) :
    0 < B u u := by
  rcases (hpsd u).lt_or_eq with h | h
  · exact h
  · exact absurd (htrivial u h.symm) hu
