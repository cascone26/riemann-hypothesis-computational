import PSDNullRadical

/-!
# §1745 — The radical of a PSD form is a linear subspace

Upgrading §1744: the null set `{u | B u u = 0}` of a positive-semidefinite symmetric bilinear form is closed
under addition and scalar multiplication — it is a submodule.  (Closure under `+` uses §1744: for two null
vectors the cross term `B a b` also vanishes, so `B (a+b) (a+b) = 0`.)

For the PSD-saturated Weil form this says the verified null space is a genuine *linear subspace* — the radical
`rad(W)` — on which the form vanishes identically.  Weil positivity is thus the statement that `W` is positive
definite on the quotient `M / rad(W)`, with the saturation `Λ=0` recorded as `rad(W) ≠ 0`.  EXIT0 / std-axioms.
-/

/-- **The radical is a submodule.**  The null set of a PSD symmetric bilinear form is a linear subspace. -/
def psdRadical {M : Type*} [AddCommGroup M] [Module ℝ M]
    (B : LinearMap.BilinForm ℝ M) (hsymm : ∀ x y, B x y = B y x) (hpsd : ∀ x, 0 ≤ B x x) :
    Submodule ℝ M where
  carrier := {u | B u u = 0}
  zero_mem' := by simp
  add_mem' := by
    intro a b ha hb
    simp only [Set.mem_setOf_eq] at ha hb ⊢
    have h1 : B a b = 0 := psd_null_in_radical B hsymm hpsd a ha b
    have h2 : B b a = 0 := by rw [hsymm b a]; exact h1
    simp only [map_add, LinearMap.add_apply, ha, hb, h1, h2, add_zero, zero_add]
  smul_mem' := by
    intro c a ha
    simp only [Set.mem_setOf_eq] at ha ⊢
    simp only [map_smul, LinearMap.smul_apply, smul_eq_mul, ha, mul_zero]
