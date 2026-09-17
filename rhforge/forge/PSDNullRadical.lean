import PSDCauchySchwarz

/-!
# §1744 — Null vectors of a PSD form are radical (degenerate) directions

Immediate from PSD Cauchy–Schwarz §1743: if the archimedean-minus-prime Weil form is positive semidefinite and a
test vector has zero self-energy, it pairs to zero with everything:

`B` symmetric, `∀ x, 0 ≤ B x x`, `B u u = 0`  ⟹  `∀ w, B u w = 0`.

This is exactly the structure tonight's numerics exhibit: the PSD-saturated Weil form (`λ_min = 0`, verified) has
a *growing low-frequency null space*, and this brick shows those null modes are genuine radical vectors — the
degenerate directions of the form, not merely small.  At the `Λ=0` knife-edge the saturation lives precisely in
this radical.  EXIT0 / std-axioms.
-/

/-- **Null ⟹ radical.**  A zero-self-energy vector of a PSD symmetric form pairs to zero with all vectors. -/
theorem psd_null_in_radical {M : Type*} [AddCommGroup M] [Module ℝ M]
    (B : LinearMap.BilinForm ℝ M) (hsymm : ∀ x y, B x y = B y x)
    (hpsd : ∀ x, 0 ≤ B x x) (u : M) (hu : B u u = 0) (w : M) :
    B u w = 0 := by
  have hcs := psd_cauchy_schwarz B hsymm hpsd u w
  rw [hu, zero_mul] at hcs
  nlinarith [sq_nonneg (B u w), hcs]
