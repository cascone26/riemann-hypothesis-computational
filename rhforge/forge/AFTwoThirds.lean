import Mathlib
import AFProportion

/-!
# Alpöge–Furman two-thirds (RH session 2026-08-29)

The unconditional `2/3` specialisation of `af_proportion` at Montgomery's second-moment value `R = 4/3`
(`R(ψ₀) = 4/3` for the indicator window): `(2/3)·N ≤ rank P₁ = #{pᵢ ≠ 0}`.
-/

open Finset

/-- At the indicator-window value `R = 4/3`: `(2/3)·N ≤ #{pᵢ ≠ 0}`. -/
theorem af_two_thirds {ι : Type*} (s : Finset ι) (p q : ι → ℝ) (N : ℝ)
    (hp : ∀ i ∈ s, 0 ≤ p i)
    (hN : (∑ i ∈ s, p i) + (∑ i ∈ s, q i) = N)
    (hR : (∑ i ∈ s, (p i + q i) ^ 2) = (4/3) * N)
    (hb : (∑ i ∈ s, p i) + 2 * ((s.filter (fun i => 0 < q i)).card : ℝ) ≤ N) :
    (2/3) * N ≤ ((s.filter (fun i => p i ≠ 0)).card : ℝ) := by
  have h := af_proportion s p q N (4/3) hp hN hR hb
  linarith
