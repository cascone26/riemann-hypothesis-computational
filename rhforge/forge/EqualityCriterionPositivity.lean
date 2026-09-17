import PSDPosDefTrivialRadical
import Mathlib.Tactic

/-!
# §1754 — Suzuki's equality criterion ⟹ Weil positivity + anisotropy (the bridge)

Closes the loop from Suzuki's RH equality-criterion (Thm 1.4/Cor 1.5: `⟨v,v⟩_W = c·I v v` on `V°(0)`, `c>0`,
`I` the L²-inner-product) to the positivity/radical framework (§1748).  If the Weil form equals a positive
multiple of a positive-semidefinite reference form on the diagonal, then:

* **positivity**: `W` is PSD (its self-energy is `c ≥ 0` times a nonnegative quantity);
* **anisotropy from strict `I`**: if additionally `I` is anisotropic (`I v v = 0 → v = 0`), then `W` is positive
  DEFINITE — no radical, RH-positivity in its strict (nondegenerate) form.

Since Suzuki's reference `I` is the genuine L²-inner-product (PD, anisotropic — exactly what §1738 establishes
for the archimedean metric on the eigenbasis), the equality criterion delivers strict positivity directly.  So
`equality on V°(0) ⟹ RH` factors through this brick + §1748.  EXIT0 / std-axioms.
-/

variable {M : Type*} [AddCommGroup M] [Module ℝ M]

/-- **Equality-criterion ⟹ PSD.**  If `W v v = c · I v v` with `c ≥ 0` and `I` PSD, then `W` is PSD. -/
theorem weil_psd_of_equality (W I : LinearMap.BilinForm ℝ M) (c : ℝ)
    (hc : 0 ≤ c) (hIpsd : ∀ v, 0 ≤ I v v) (hdiag : ∀ v, W v v = c * I v v) (v : M) :
    0 ≤ W v v := by
  rw [hdiag v]; exact mul_nonneg hc (hIpsd v)

/-- **Equality-criterion ⟹ positive definite.**  If moreover `c > 0` and `I` is anisotropic, `W` is PD. -/
theorem weil_posdef_of_equality (W I : LinearMap.BilinForm ℝ M) (c : ℝ)
    (hc : 0 < c) (hIpsd : ∀ v, 0 ≤ I v v) (hIaniso : ∀ v, I v v = 0 → v = 0)
    (hdiag : ∀ v, W v v = c * I v v) (v : M) (hv : v ≠ 0) :
    0 < W v v := by
  rw [hdiag v]
  have : 0 < I v v := lt_of_le_of_ne (hIpsd v) (fun h => hv (hIaniso v h.symm))
  exact mul_pos hc this
