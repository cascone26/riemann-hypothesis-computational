import NormFactorLe
import Mathlib.Analysis.Complex.AbsMax

set_option maxHeartbeats 1000000

open Complex Metric

/-!
# §1673 ★★ — JENSEN-LITE ZERO COUNT: `3^(Σ mult) ≤ M/‖f c‖` — zeros in the quarter-ball

The zero-counting input of the Landau ladder, by the elementary radius-gap argument (no integrals):
if `f = P·g` with `P` the polynomial of zeros in `closedBall c (R/4)` (multiplicities `m`), `f`
bounded by `M` on `sphere c R`, and `g` analytic, then each zero contributes a factor `3`:
on the sphere `‖s-ρ‖ ≥ 3R/4` so `‖g‖ ≤ M/(3R/4)^Σm` there, max-modulus pushes that to `c`, while
`‖P c‖ ≤ (R/4)^Σm`; multiplying, `‖f c‖·3^Σm ≤ M`. So `Σ m ≤ log(M/‖f c‖)/log 3` — the count is
`O(log M)`. Feeds §1671's `(Σ|d|)/δ` term.
-/

theorem three_pow_sum_mult_le_of_factor {c : ℂ} {R M : ℝ} {f g P : ℂ → ℂ} {S : Finset ℂ}
    {m : ℂ → ℕ} (hR : 0 < R) (hPdef : P = fun s => ∏ ρ ∈ S, (s - ρ) ^ (m ρ))
    (hS : ∀ ρ ∈ S, ρ ∈ closedBall c (R / 4))
    (hfac : ∀ s ∈ sphere c R, f s = P s * g s) (hfacc : f c = P c * g c)
    (hfsphere : ∀ s ∈ sphere c R, ‖f s‖ ≤ M)
    (hgdiff : DiffContOnCl ℂ g (ball c R)) (hfc0 : f c ≠ 0) :
    (3 : ℝ) ^ (∑ ρ ∈ S, m ρ) ≤ M / ‖f c‖ := by
  set n : ℕ := ∑ ρ ∈ S, m ρ with hndef
  have hfcpos : 0 < ‖f c‖ := norm_pos_iff.mpr hfc0
  have h34 : (0 : ℝ) < 3 * R / 4 := by linarith
  -- lower bound for ‖P‖ on the sphere: each factor ≥ 3R/4
  have hPlow : ∀ s ∈ sphere c R, (3 * R / 4) ^ n ≤ ‖P s‖ := by
    intro s hs
    rw [hPdef, norm_prod]
    simp_rw [norm_pow]
    rw [hndef, ← Finset.prod_pow_eq_pow_sum]
    apply Finset.prod_le_prod
    · intro ρ _
      positivity
    · intro ρ hρ
      have hbase : 3 * R / 4 ≤ ‖s - ρ‖ := by
        have hρc : dist ρ c ≤ R / 4 := mem_closedBall.mp (hS ρ hρ)
        have hsc : dist s c = R := mem_sphere.mp hs
        have htri := dist_triangle s ρ c
        rw [← dist_eq_norm]
        linarith
      gcongr
  -- sphere bound for g, then max modulus to the center
  have hgsphere : ∀ s ∈ sphere c R, ‖g s‖ ≤ M / (3 * R / 4) ^ n := by
    intro s hs
    exact norm_factor_le (hfac s hs) (by positivity) (hPlow s hs) (hfsphere s hs)
  have hgc : ‖g c‖ ≤ M / (3 * R / 4) ^ n := by
    have hfr : ∀ u ∈ frontier (ball c R), ‖g u‖ ≤ M / (3 * R / 4) ^ n := by
      intro u hu
      rw [frontier_ball c hR.ne'] at hu
      exact hgsphere u hu
    exact Complex.norm_le_of_forall_mem_frontier_norm_le isBounded_ball hgdiff hfr
      (subset_closure (mem_ball_self hR))
  -- upper bound for ‖P c‖: each factor ≤ R/4
  have hPc : ‖P c‖ ≤ (R / 4) ^ n := by
    rw [hPdef, norm_prod]
    simp_rw [norm_pow]
    rw [hndef, ← Finset.prod_pow_eq_pow_sum]
    apply Finset.prod_le_prod
    · intro ρ _
      positivity
    · intro ρ hρ
      have hbase : ‖c - ρ‖ ≤ R / 4 := by
        rw [← dist_eq_norm, dist_comm c ρ]
        exact mem_closedBall.mp (hS ρ hρ)
      gcongr
  -- combine, cross-multiplied: ‖f c‖·3^n ≤ M
  have hXpos : (0 : ℝ) < (3 * R / 4) ^ n := by positivity
  have hkey : ‖f c‖ * 3 ^ n ≤ M := by
    calc ‖f c‖ * 3 ^ n = ‖P c‖ * ‖g c‖ * 3 ^ n := by rw [hfacc, norm_mul]
      _ ≤ (R / 4) ^ n * (M / (3 * R / 4) ^ n) * 3 ^ n := by
          gcongr
      _ = (R / 4 * 3) ^ n * (M / (3 * R / 4) ^ n) := by rw [mul_pow]; ring
      _ = (3 * R / 4) ^ n * (M / (3 * R / 4) ^ n) := by
          rw [show R / 4 * 3 = 3 * R / 4 by ring]
      _ = M / (3 * R / 4) ^ n * (3 * R / 4) ^ n := by ring
      _ = M := div_mul_cancel₀ M hXpos.ne'
  rw [le_div_iff₀ hfcpos, mul_comm]
  exact hkey

#print axioms three_pow_sum_mult_le_of_factor
