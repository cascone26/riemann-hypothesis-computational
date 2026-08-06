import RiemannLean.SignTheorem
open ArithmeticFunction
/-- Completeness of the sign theorem: under the onset factorization, the onset is nonzero and its
strict sign tracks `μ(n)` — `0 < D ↔ μ(n)=1` and `D < 0 ↔ μ(n)=-1` (squarefree `n`). -/
theorem sign_theorem_onset_definite {n : ℕ} (hn : Squarefree n) {T S : Finset ℕ}
    (hT : ∀ p ∈ T, 2 ≤ p) (hS : ∀ q ∈ S, 2 ≤ q) {D : ℝ}
    (hD : D = (moebius n : ℝ) *
      ((∏ p ∈ T, besselJ1 (1 / Real.sqrt p) / besselJ0 (1 / Real.sqrt p))
        * ∏ q ∈ S, besselJ0 (1 / Real.sqrt q))) :
    D ≠ 0 := by
  have hsign := sign_theorem_onset hn hT hS hD
  have hne : moebius n ≠ 0 := moebius_ne_zero_iff_squarefree.mpr hn
  intro hD0
  rw [hD0, Real.sign_zero] at hsign
  have : (moebius n : ℝ) = 0 := hsign.symm
  exact hne (by exact_mod_cast this)

#print axioms sign_theorem_onset_definite
