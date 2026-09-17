import Mathlib
import RankTraceCodiag

/-!
# Alpöge–Furman proportion certificate (RH session 2026-08-29)

The AF headline `rank P₁ ≥ (2 − R)·N` (eq chain lines 133–137), derived in-forge from the co-diagonal
rank–trace inequality plus the three analytic inputs as hypotheses:
  `∑(pᵢ+qᵢ) = N`   (trace normalisation `tr Γ̃ = N`),
  `∑(pᵢ+qᵢ)² = R·N` (`‖Γ̃‖²_HS = R·N`, Montgomery's 2nd moment),
  `∑pᵢ + 2·#{qᵢ>0} ≤ N` (AF's counting `tr P₁ + 2·npl(Q') ≤ N`).
Conclusion: `(2 − R)·N ≤ rank P₁ = #{pᵢ≠0}`.  At `R = 4/3` this is the unconditional `2/3`.
-/

open Finset

/-- AF proportion certificate: `(2 − R)·N ≤ #{pᵢ ≠ 0}`. -/
theorem af_proportion {ι : Type*} (s : Finset ι) (p q : ι → ℝ) (N R : ℝ)
    (hp : ∀ i ∈ s, 0 ≤ p i)
    (hN : (∑ i ∈ s, p i) + (∑ i ∈ s, q i) = N)
    (hR : (∑ i ∈ s, (p i + q i) ^ 2) = R * N)
    (hb : (∑ i ∈ s, p i) + 2 * ((s.filter (fun i => 0 < q i)).card : ℝ) ≤ N) :
    (2 - R) * N ≤ ((s.filter (fun i => p i ≠ 0)).card : ℝ) := by
  have key := rank_trace_codiagonal s p q hp
  rw [hR] at key
  nlinarith [key, hN, hb]
