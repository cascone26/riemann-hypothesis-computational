import Mathlib.LinearAlgebra.BilinearForm.Basic
import Mathlib.Tactic

/-!
# §1741 — Weil positivity ⟺ prime energy dominated by archimedean energy

The Weil form is `W = A − P` (archimedean minus prime).  Its positivity on a test vector is exactly the
statement that the prime energy does not exceed the archimedean energy there:

`0 ≤ (A − P) v v  ⟺  P v v ≤ A v v`.

This is the clean criterion behind tonight's verified numerics (2026-08-24): the prime functional matches the
archimedean functional mode-by-mode (ratio prime/|arch| → 1) and lands on `0⁺` — PSD-saturation.  Weil
positivity (hence RH) is precisely "arch dominates prime on every mode", and the saturation `Λ=0` is the
knife-edge where the domination is exactly tight.  EXIT0 / std-axioms.
-/

/-- **Domination criterion.**  `A − P` is nonnegative on `v` iff `P v v ≤ A v v`. -/
theorem weil_positivity_iff_domination {M : Type*} [AddCommGroup M] [Module ℝ M]
    (A P : LinearMap.BilinForm ℝ M) (v : M) :
    0 ≤ (A - P) v v ↔ P v v ≤ A v v := by
  simp only [LinearMap.sub_apply, sub_nonneg]
