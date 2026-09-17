import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.Analysis.Meromorphic.Divisor

/-!
# §1415 — centre product upper bound `‖∏ (c-ρ)^{m ρ}‖ ≤ r^{Σ m ρ}`

Mirror of §1375/§1376 (which lower-bounded `‖P_half‖` on the sphere): at the **centre** `c`, every half-disk zero
`ρ ∈ closedBall c r` satisfies `‖c-ρ‖ ≤ r`, so the weighted product is `≤ r^(Σ m ρ)`.  Combined with the
centre factorization `ζ c = P_half(c)·g_full(c)`, this gives the LOWER bound `‖g_full c‖ ≥ ‖ζ c‖/r^N` — the missing
half that makes the `N`-terms in §1402's `BOUND` cancel, leaving `BOUND ≤ 4·log(Mζ/‖ζ c‖)/(R/2)` (explicit, no
`g₀`).  EXIT0 / std-axioms.
-/

theorem norm_prod_sub_center_le {S : Finset ℂ} {c : ℂ} {r : ℝ} (hr : 0 ≤ r) (m : ℂ → ℕ)
    (hS : ∀ ρ ∈ S, ρ ∈ Metric.closedBall c r) :
    ‖∏ ρ ∈ S, (c - ρ) ^ (m ρ)‖ ≤ r ^ (∑ ρ ∈ S, m ρ) := by
  rw [norm_prod, ← Finset.prod_pow_eq_pow_sum]
  apply Finset.prod_le_prod
  · intro ρ _; exact norm_nonneg _
  · intro ρ hρ
    rw [norm_pow]
    apply pow_le_pow_left₀ (norm_nonneg _) _
    have h := hS ρ hρ
    rw [Metric.mem_closedBall, Complex.dist_eq] at h
    rwa [norm_sub_rev]

#print axioms norm_prod_sub_center_le
