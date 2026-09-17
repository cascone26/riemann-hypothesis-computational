import Mathlib.Analysis.SpecialFunctions.Log.Basic

/-!
# §1434 — the `hfar` threshold: `σK - R/2 ≤ 1 - 1/(14K)` for `8 ≤ 7·K·R`

Pure scalar lemma discharging §1433's `hfar`.  With `σK = 1 + 1/(2K)`, the inequality reduces to
`1/(2K) + 1/(14K) = 8/(14K) ≤ R/2`, i.e. `8 ≤ 7·K·R`.  For `R = 1/4` this is `K ≥ 32/7 ≈ 4.58` — always satisfied
in the fixed point since `K = 3A+4B+C ≍ log|γ| → ∞`.  EXIT0 / std-axioms.
-/

theorem hfar_holds {K R : ℝ} (hK : 0 < K) (h : 8 ≤ 7 * K * R) :
    (1 + 1 / (2 * K)) - R / 2 ≤ 1 - 1 / (14 * K) := by
  have h14K : (0 : ℝ) < 14 * K := by linarith
  have key : 1 / (2 * K) + 1 / (14 * K) ≤ R / 2 := by
    have heq : 1 / (2 * K) + 1 / (14 * K) = 8 / (14 * K) := by field_simp; ring
    rw [heq, div_le_iff₀ h14K]
    nlinarith [h]
  linarith [key]

#print axioms hfar_holds
