import Mathlib.LinearAlgebra.BilinearForm.Basic
import Mathlib.Tactic

/-!
# §1758 — Birman–Schwinger gap bound for the Schrödinger reframing of λa

The Dirichlet-form reframing writes the localized Weil form as `A − P` (positive kinetic `A` minus potential
`P`), so `λa` is a Schrödinger ground-state energy.  The quantitative Birman–Schwinger / relative-boundedness
criterion: if the potential is *sub-unit-coupled* to the kinetic energy — `P v v ≤ (1−ε)·A v v` for all `v`,
some `ε ≥ 0` — then the Schrödinger form has a proportional positive gap:

`P v v ≤ (1−ε)·A v v` and `A` PSD  ⟹  `ε·A v v ≤ (A − P) v v`.

In particular `ε > 0` with `A` positive definite gives `A − P` positive definite (`λa > 0`).  This is the exact
tool for the frontier `RH ⟺ λa > 0 ∀a`: bound the prime-potential coupling to the archimedean kinetic energy
strictly below `1`, uniformly in `a`.  EXIT0 / std-axioms.
-/

variable {M : Type*} [AddCommGroup M] [Module ℝ M]

/-- **Birman–Schwinger gap.**  Sub-unit coupling of the potential gives a proportional positive gap. -/
theorem birman_schwinger_gap (A P : LinearMap.BilinForm ℝ M) (ε : ℝ)
    (hApsd : ∀ v, 0 ≤ A v v) (hcoup : ∀ v, P v v ≤ (1 - ε) * A v v) (v : M) :
    ε * A v v ≤ (A - P) v v := by
  have h := hcoup v
  have : (A - P) v v = A v v - P v v := by simp [LinearMap.sub_apply]
  rw [this]; nlinarith [h]

/-- **Strict positivity from sub-unit coupling.**  `ε > 0`, `A` PD, sub-unit coupling ⟹ `A − P` is PD. -/
theorem schrodinger_posdef_of_subunit (A P : LinearMap.BilinForm ℝ M) (ε : ℝ)
    (hε : 0 < ε) (hApd : ∀ v, v ≠ 0 → 0 < A v v) (hApsd : ∀ v, 0 ≤ A v v)
    (hcoup : ∀ v, P v v ≤ (1 - ε) * A v v) (v : M) (hv : v ≠ 0) :
    0 < (A - P) v v :=
  lt_of_lt_of_le (mul_pos hε (hApd v hv)) (birman_schwinger_gap A P ε hApsd hcoup v)
