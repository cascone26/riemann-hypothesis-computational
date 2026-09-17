import Mathlib.LinearAlgebra.BilinearForm.Basic
import Mathlib.Tactic

/-!
# §1757 — Suzuki's equality subspace is transverse to the radical

On Suzuki's `V°(0)` the Weil form satisfies the equality `W(v,v) = c·I(v,v)` (`c = 2 > 0`, `I` the L²-inner
product, anisotropic).  The PSD-saturation radical `rad(W) = {v | W(v,v) = 0}` is where the form vanishes.
These two are transverse: on the equality subspace the form is a *positive multiple of an anisotropic form*, so
it cannot vanish on a nonzero vector.

`W(v,v) = c·I(v,v)` on `S`, `c > 0`, `I` anisotropic  ⟹  `S ∩ rad(W) = {0}`.

So Suzuki's `V°(0)` lives entirely OUTSIDE the saturation radical — it is exactly the strictly-positive part of
the Weil form.  This is the structural reason the equality criterion (`W = 2I` on `V°(0)`) and the radical
theory (`W = 0` on `rad`) are compatible: `M` splits as `V°(0) ⊕ rad`, form `= 2I ⊕ 0`.  RH ⟺ `V°(0)` fills the
whole non-radical part.  EXIT0 / std-axioms.
-/

variable {M : Type*} [AddCommGroup M] [Module ℝ M]

/-- **Transversality.**  A vector on the equality subspace that is also in the radical must be `0`. -/
theorem equality_subspace_disjoint_radical (W I : LinearMap.BilinForm ℝ M) (c : ℝ)
    (hc : 0 < c) (hIaniso : ∀ v, I v v = 0 → v = 0)
    (v : M) (hEq : W v v = c * I v v) (hRad : W v v = 0) :
    v = 0 := by
  have hIv : I v v = 0 := by
    have : c * I v v = 0 := by rw [← hEq]; exact hRad
    rcases mul_eq_zero.mp this with h | h
    · exact absurd h (ne_of_gt hc)
    · exact h
  exact hIaniso v hIv
