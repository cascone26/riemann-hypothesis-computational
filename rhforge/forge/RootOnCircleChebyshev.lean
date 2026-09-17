import JoukowskiOnCircle
import Mathlib.RingTheory.Polynomial.Chebyshev

/-!
# §1785 — Capstone: a palindromic root is on the unit circle ⟺ its Chebyshev image is real in [-2,2]

Assembles the Chebyshev on-circle route (§1781 Joukowski criterion + §1784 degree-reduction identity
`P(ρ) = ρᵐ·Q(ρ+ρ⁻¹)`).  For a nonzero `ρ`:

  `ρ` is a root of `P` lying ON the unit circle  ⟺  `s = ρ+ρ⁻¹` is a REAL root of `Q` in `[-2,2]`.

This closes the on-circle criterion by REAL-ROOTEDNESS of the half-degree polynomial `Q` (Sturm-tractable),
discharging §1777's owed input ("self-inversive PSD-singular ⟹ zeros on `|z|=1`") without the mathlib-absent
Herglotz / trigonometric-moment theory.  EXIT0 / std-axioms.
-/

open Complex Polynomial

/-- **On-circle root criterion via Chebyshev reduction.**  Given the degree-reduction factorization
`P.eval ρ = ρᵐ · Q.eval (ρ+ρ⁻¹)` (§1784) and `ρ ≠ 0`: `ρ` is a root of `P` on the unit circle iff
`ρ+ρ⁻¹` is a real root of `Q` with `|ρ+ρ⁻¹| ≤ 2`. -/
theorem root_on_circle_iff_chebyshev {P Q : ℂ[X]} {ρ : ℂ} (m : ℕ) (hρ : ρ ≠ 0)
    (hfact : P.eval ρ = ρ ^ m * Q.eval (ρ + ρ⁻¹)) :
    (P.eval ρ = 0 ∧ normSq ρ = 1) ↔
      (Q.eval (ρ + ρ⁻¹) = 0 ∧ (ρ + ρ⁻¹).im = 0 ∧ normSq (ρ + ρ⁻¹) ≤ 4) := by
  have hpm : ρ ^ m ≠ 0 := pow_ne_zero m hρ
  have hP : P.eval ρ = 0 ↔ Q.eval (ρ + ρ⁻¹) = 0 := by
    rw [hfact, mul_eq_zero]
    exact ⟨fun h => h.resolve_left hpm, fun h => Or.inr h⟩
  have hcirc := joukowski_on_circle ρ hρ
  rw [hP]
  constructor
  · rintro ⟨hq, hns⟩; exact ⟨hq, hcirc.mpr hns⟩
  · rintro ⟨hq, hc⟩; exact ⟨hq, hcirc.mp hc⟩
