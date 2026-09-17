import Mathlib.LinearAlgebra.BilinearForm.Basic
import Mathlib.Tactic

/-!
# §1753 — Polarization: a diagonal form-equals-norm criterion IS a full isometry

Suzuki's unconditional RH criterion (Thm 1.4 / Cor 1.5, arXiv:2301.00421) is stated as a DIAGONAL equality:
`⟨ψ, ψ⟩_W = 2‖ψ‖²` for all `ψ` in a subspace `V°(0)`.  Its deep content is the full ISOMETRY
`⟨u, v⟩_W = 2⟨u, v⟩` (the de Branges structure) — and polarization shows the diagonal equality already
CONTAINS the bilinear one:

for symmetric bilinear forms `W`, `I` and a scalar `c`, if `W v v = c · I v v` for all `v`, then
`W u v = c · I u v` for all `u, v`.

So Suzuki's checkable diagonal criterion is equivalent to the Weil form being `c`-times an inner product on
`V°(0)` — exactly the "Weil form ≅ L²" isometry that makes `HW` a Hilbert–Pólya / de Branges space.  This is
the abstract bridge from the equality-criterion (diagonal, §1735–1738 shape) to the isometry (bilinear).
EXIT0 / std-axioms.
-/

variable {M : Type*} [AddCommGroup M] [Module ℝ M]

/-- **Polarization of a form-equals-norm identity.**  If two symmetric bilinear forms agree on the diagonal
up to a scalar `c`, they agree as bilinear forms. -/
theorem bilin_eq_of_diag_eq (W I : LinearMap.BilinForm ℝ M) (c : ℝ)
    (hWsymm : ∀ x y, W x y = W y x) (hIsymm : ∀ x y, I x y = I y x)
    (hdiag : ∀ v, W v v = c * I v v) (u v : M) :
    W u v = c * I u v := by
  have hW : W (u + v) (u + v) = W u u + 2 * W u v + W v v := by
    simp only [map_add, LinearMap.add_apply]; rw [hWsymm v u]; ring
  have hI : I (u + v) (u + v) = I u u + 2 * I u v + I v v := by
    simp only [map_add, LinearMap.add_apply]; rw [hIsymm v u]; ring
  have huv := hdiag (u + v)
  rw [hW, hI] at huv
  have hu := hdiag u
  have hv := hdiag v
  nlinarith [huv, hu, hv]
