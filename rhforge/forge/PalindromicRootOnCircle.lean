import PalindromicDicksonFactor
import JoukowskiOnCircle

/-!
# §1789 — End-to-end on-circle criterion for palindromic polynomials (headline theorem)

Composes the whole Chebyshev chain (§1788 factorization + §1781 Joukowski criterion) into a single
statement.  For a palindromic coefficient function `v` (`v i = v (2m-i)`) and `ρ ≠ 0`:

  `ρ` is a root of the degree-`2m` palindromic polynomial lying ON the unit circle
    ⟺
  `s = ρ + ρ⁻¹` is a REAL root (in `[-2,2]`) of the half-degree Dickson combination `Q(s)`.

This is the complete, unconditional, mathlib-absent reduction of "palindromic-polynomial zeros on `|z|=1`"
to "real-rootedness of a half-degree polynomial in `[-2,2]`", verifying that the 8-brick chain §1781–1788
composes end to end.  EXIT0 / std-axioms.
-/

open Complex Polynomial

/-- **Palindromic on-circle criterion (headline).**  A root `ρ ≠ 0` of a palindromic polynomial lies on the
unit circle iff its Chebyshev image `ρ+ρ⁻¹` is a real root of the Dickson combination `Q` with `|ρ+ρ⁻¹| ≤ 2`. -/
theorem palindromic_root_on_circle_iff (v : ℕ → ℂ) (ρ : ℂ) (hρ : ρ ≠ 0) (m : ℕ)
    (hpal : ∀ i, i ≤ 2 * m → v i = v (2 * m - i)) :
    ((∑ i ∈ Finset.range (2 * m + 1), v i * ρ ^ i) = 0 ∧ normSq ρ = 1) ↔
      ((v m + ∑ k ∈ Finset.range m, v (m + 1 + k) * (dickson (k + 1)).eval (ρ + ρ⁻¹)) = 0
        ∧ (ρ + ρ⁻¹).im = 0 ∧ normSq (ρ + ρ⁻¹) ≤ 4) := by
  have hfact := palindromic_dickson_factor v ρ hρ m hpal
  set Q := v m + ∑ k ∈ Finset.range m, v (m + 1 + k) * (dickson (k + 1)).eval (ρ + ρ⁻¹) with hQ
  have hpm : ρ ^ m ≠ 0 := pow_ne_zero m hρ
  have hsum0 : ((∑ i ∈ Finset.range (2 * m + 1), v i * ρ ^ i) = 0) ↔ Q = 0 := by
    rw [hfact, mul_eq_zero]
    exact ⟨fun h => h.resolve_left hpm, fun h => Or.inr h⟩
  have hcirc := joukowski_on_circle ρ hρ
  rw [hsum0]
  constructor
  · rintro ⟨hq, hns⟩; exact ⟨hq, hcirc.mpr hns⟩
  · rintro ⟨hq, hc⟩; exact ⟨hq, hcirc.mp hc⟩
