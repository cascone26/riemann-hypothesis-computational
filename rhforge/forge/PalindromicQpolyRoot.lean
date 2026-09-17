import PalindromicDicksonFactor
import QpolyDef

/-!
# §1814 — Nonzero root of P ⟺ Chebyshev image is a root of `Qpoly`

Combines §1788 (factorization `P(ρ) = ρᵐ·Q(ρ+ρ⁻¹)`) with §1813 (`Qpoly`) to give the root correspondence:
for `ρ ≠ 0`, `ρ` is a root of the palindromic `P` iff `ρ+ρ⁻¹` is a root of the half-degree `Qpoly`.  This
maps P's nonzero roots onto `Qpoly`'s roots via the Joukowski map, the bridge that lets `(Qpoly v m).roots`
(a `Multiset ℂ`) carry the on-circle information into the real-rootedness iffs §1806/§1811.  EXIT0 / std-axioms.
-/

open Polynomial

/-- **Root correspondence.**  `ρ≠0` is a root of `P` iff `ρ+ρ⁻¹` is a root of `Qpoly`. -/
theorem palindromic_root_iff_Qpoly_root (v : ℕ → ℂ) (ρ : ℂ) (hρ : ρ ≠ 0) (m : ℕ)
    (hpal : ∀ i, i ≤ 2 * m → v i = v (2 * m - i)) :
    (∑ i ∈ Finset.range (2 * m + 1), v i * ρ ^ i) = 0 ↔ (Qpoly v m).IsRoot (ρ + ρ⁻¹) := by
  rw [palindromic_dickson_factor v ρ hρ m hpal, ← Qpoly_eval]
  simp only [Polynomial.IsRoot.def]
  constructor
  · intro h; exact (mul_eq_zero.mp h).resolve_left (pow_ne_zero m hρ)
  · intro h; rw [h, mul_zero]
