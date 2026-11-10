import Mathlib

/-!
# S3037 — radius gap avoiding a finite set (lane 1, Poisson–Landau brick P2a)

`exists_gap`: for any finite T ⊆ ℝ and a < b there are a < r < r' < b with every t ∈ T outside [r, r'].
Use: T = {‖u − c‖ : u zero of ζ in the big ball}; then no zero lies in the closed annulus r ≤ |z − c| ≤ r', so after
dividing out the Blaschke product of the zeros in ball c r, G = ζ/B is zero-free on ball c r' (r < r'), analytic there,
and |G| = |ζ| on sphere c r — exactly the hypotheses of the Schwarz formula (S3034) at radius r.  The leading Poisson–Landau
coefficient (2/r)(2π)⁻¹(r/2)(π/2) = 1/4 does not depend on r, so any r ∈ (a, b) works.
-/

namespace ZeroGap

theorem exists_gap (T : Finset ℝ) {a b : ℝ} (hab : a < b) :
    ∃ r r' : ℝ, a < r ∧ r < r' ∧ r' < b ∧ ∀ t ∈ T, t < r ∨ r' < t := by
  induction T using Finset.induction_on with
  | empty =>
    exact ⟨a + (b - a) / 3, a + 2 * (b - a) / 3, by linarith, by linarith, by linarith, by simp⟩
  | insert t T ht ih =>
    obtain ⟨r, r', har, hrr, hrb, hT⟩ := ih
    by_cases hout : t < r ∨ r' < t
    · refine ⟨r, r', har, hrr, hrb, ?_⟩
      intro x hx
      rcases Finset.mem_insert.mp hx with rfl | hx
      · exact hout
      · exact hT x hx
    · push Not at hout
      obtain ⟨h1, h2⟩ := hout
      by_cases hmid : (r' - r) / 2 ≤ t - r
      · -- gap (r, t) is long enough
        have htr : 0 < t - r := by linarith
        refine ⟨r + (t - r) / 3, r + 2 * (t - r) / 3, by linarith, by linarith, by linarith, ?_⟩
        intro x hx
        rcases Finset.mem_insert.mp hx with rfl | hx
        · right; linarith
        · rcases hT x hx with h | h
          · left; linarith
          · right; linarith
      · push Not at hmid
        have htr : 0 < r' - t := by linarith
        refine ⟨t + (r' - t) / 3, t + 2 * (r' - t) / 3, by linarith, by linarith, by linarith, ?_⟩
        intro x hx
        rcases Finset.mem_insert.mp hx with rfl | hx
        · left; linarith
        · rcases hT x hx with h | h
          · left; linarith
          · right; linarith

end ZeroGap

#print axioms ZeroGap.exists_gap
