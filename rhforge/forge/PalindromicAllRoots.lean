import PalindromicRootOnCircle
import PalindromicDicksonFactor

/-!
# §1812 — All-roots on-circle criterion (lift of §1789 over the root set)

Lifts the per-root iff §1789 to a statement over ALL nonzero roots of the palindromic polynomial:
every root lies on the unit circle IFF every root's Chebyshev image `ρ+ρ⁻¹` is real and in `[−2,2]`.
Uses §1788 (`palindromic_dickson_factor`, so `P(ρ)=0 ⟺ Q(ρ+ρ⁻¹)=0` since `ρᵐ≠0`) to feed §1789 pointwise.
This is the root-set form that connects the on-circle module to the real-rootedness iffs (§1806/§1811)
applied to `Q`.  EXIT0 / std-axioms.
-/

open Complex Polynomial

/-- **All-roots on-circle criterion.** -/
theorem palindromic_all_roots_on_circle (v : ℕ → ℂ) (m : ℕ)
    (hpal : ∀ i, i ≤ 2 * m → v i = v (2 * m - i)) :
    (∀ ρ : ℂ, ρ ≠ 0 → (∑ i ∈ Finset.range (2 * m + 1), v i * ρ ^ i) = 0 → normSq ρ = 1)
      ↔ (∀ ρ : ℂ, ρ ≠ 0 → (∑ i ∈ Finset.range (2 * m + 1), v i * ρ ^ i) = 0 →
          ((ρ + ρ⁻¹).im = 0 ∧ normSq (ρ + ρ⁻¹) ≤ 4)) := by
  constructor
  · intro h ρ hρ hroot
    exact ((palindromic_root_on_circle_iff v ρ hρ m hpal).mp ⟨hroot, h ρ hρ hroot⟩).2
  · intro h ρ hρ hroot
    have h1 := h ρ hρ hroot
    have hfact := palindromic_dickson_factor v ρ hρ m hpal
    have hQ0 : (v m + ∑ k ∈ Finset.range m, v (m + 1 + k) * (dickson (k + 1)).eval (ρ + ρ⁻¹)) = 0 := by
      rw [hfact] at hroot
      exact (mul_eq_zero.mp hroot).resolve_left (pow_ne_zero m hρ)
    exact ((palindromic_root_on_circle_iff v ρ hρ m hpal).mpr ⟨hQ0, h1.1, h1.2⟩).2
