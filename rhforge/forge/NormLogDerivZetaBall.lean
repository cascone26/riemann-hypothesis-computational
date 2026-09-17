import ZetaQuarterBallFactorization
import NormLogDerivZetaMasterV2
import NormGCenterLe

set_option maxHeartbeats 2000000

open Complex Metric

/-!
# §1683 ★★★★ — THE ζ'/ζ BALL BOUND, self-contained:
`‖ζ'/ζ(z)‖ ≤ log((M₀+1)/‖ζ(c)‖) · (1/(δ·log 3) + 32/R)`

THE fully-assembled Landau bound with NO internal factorization data in the statement: on any
pole-free ball (`1 ∉ ball c R'`, `1 < Re c`), for `z` in the `R/16`-ball keeping distance `δ` from
every ζ-zero of the quarter-ball, with `‖ζ‖ ≤ M₀` on the `R`-sphere. All the §1681 factorization
data (`S`, `m`, `g`) is instantiated and then ELIMINATED: the zero count via §1674 telescopes into
`log(M/‖ζ c‖)/log 3`, and the `‖g c‖`-term telescopes via `‖g c‖ ≥ ‖ζ c‖/(R/4)^Σ` (each divided
zero pays `log 3` BACK), leaving a bound in observable quantities only. This is `ζ'/ζ = O(log²T)`:
at `c = 1+1/log T+iT`, `R = 1/4`, `M₀ = O(T)`, `‖ζ c‖ ≥ 1/O(log T)` (§1678),
`δ = c₀/(2 log(T+2))` (§1667) — the connector-integral input for effective PNT.
-/

theorem norm_logDeriv_riemannZeta_le_of_ball {c z : ℂ} {R R' M₀ δ : ℝ}
    (hR : 0 < R) (hRR' : R < R') (hδ : 0 < δ) (hc : 1 < c.re)
    (hball : (1 : ℂ) ∉ ball c R')
    (hζsphere : ∀ s ∈ sphere c R, ‖riemannZeta s‖ ≤ M₀)
    (hz : z ∈ closedBall c (R / 16))
    (hdist : ∀ ρ ∈ closedBall c (R / 4), riemannZeta ρ = 0 → δ ≤ ‖z - ρ‖) :
    ‖logDeriv riemannZeta z‖
      ≤ Real.log ((M₀ + 1) / ‖riemannZeta c‖) * (1 / (δ * Real.log 3) + 32 / R) := by
  obtain ⟨S, m, g, hSloc, hSzero, hfac, hgdiff, hgnz⟩ :=
    exists_zeta_quarterBall_factorization hR hRR' hc hball
  set n : ℕ := ∑ ρ ∈ S, m ρ with hndef
  have hζc0 : riemannZeta c ≠ 0 := riemannZeta_ne_zero_of_one_lt_re hc
  have hζcpos : 0 < ‖riemannZeta c‖ := norm_pos_iff.mpr hζc0
  have h34 : (0 : ℝ) < 3 * R / 4 := by linarith
  have hM₀0 : 0 ≤ M₀ := by
    have hpt : c + (R : ℂ) ∈ sphere c R := by
      rw [mem_sphere, dist_eq_norm]
      simp [abs_of_pos hR]
    exact le_trans (norm_nonneg _) (hζsphere _ hpt)
  have hMpos : (0 : ℝ) < M₀ + 1 := by linarith
  -- strict hgc_lt with the +1 slack
  have hgc_le := norm_g_center_le_of_quarterBall_factor hR rfl hSloc
    (fun s hs => hfac s (sphere_subset_closedBall hs)) hζsphere hgdiff
  have hgc_lt : ‖g c‖ < (M₀ + 1) / (3 * R / 4) ^ n := by
    calc ‖g c‖ ≤ M₀ / (3 * R / 4) ^ n := hgc_le
      _ < (M₀ + 1) / (3 * R / 4) ^ n := by gcongr; linarith
  -- master bound
  have hmaster := norm_logDeriv_riemannZeta_le_master' hR hδ rfl hSloc hfac
    (fun s hs => le_trans (hζsphere s hs) (by linarith)) hgdiff hgnz hgc_lt hζc0 hz
    (fun ρ hρ => hdist ρ (hSloc ρ hρ) (hSzero ρ hρ))
  -- eliminate ‖g c‖: lower bound via the center factorization
  have hgcpos : 0 < ‖g c‖ := norm_pos_iff.mpr (hgnz c (mem_ball_self (by linarith)))
  have hR4pos : (0 : ℝ) < (R / 4) ^ n := by positivity
  have hPc : ‖∏ ρ ∈ S, (c - ρ) ^ (m ρ)‖ ≤ (R / 4) ^ n := by
    rw [norm_prod]
    simp_rw [norm_pow]
    rw [hndef, ← Finset.prod_pow_eq_pow_sum]
    apply Finset.prod_le_prod
    · intro ρ _
      positivity
    · intro ρ hρ
      have hbase : ‖c - ρ‖ ≤ R / 4 := by
        rw [← dist_eq_norm, dist_comm c ρ]
        exact mem_closedBall.mp (hSloc ρ hρ)
      gcongr
  have hgc_lower : ‖riemannZeta c‖ / (R / 4) ^ n ≤ ‖g c‖ := by
    have hfacc := hfac c (mem_closedBall_self (by linarith))
    rw [div_le_iff₀ hR4pos]
    calc ‖riemannZeta c‖ = ‖∏ ρ ∈ S, (c - ρ) ^ (m ρ)‖ * ‖g c‖ := by
          rw [hfacc, norm_mul]
      _ ≤ (R / 4) ^ n * ‖g c‖ := by gcongr
      _ = ‖g c‖ * (R / 4) ^ n := by ring
  -- the log telescoping
  set L : ℝ := Real.log ((M₀ + 1) / ‖riemannZeta c‖) with hLdef
  have hterm2 : Real.log ((M₀ + 1) / (3 * R / 4) ^ n) - Real.log ‖g c‖ ≤ L := by
    have h1 : Real.log ((M₀ + 1) / (3 * R / 4) ^ n)
        = Real.log (M₀ + 1) - n * Real.log (3 * R / 4) := by
      rw [Real.log_div hMpos.ne' (by positivity), Real.log_pow]
    have h2 : Real.log (‖riemannZeta c‖ / (R / 4) ^ n) ≤ Real.log ‖g c‖ :=
      Real.log_le_log (by positivity) hgc_lower
    have h3 : Real.log (‖riemannZeta c‖ / (R / 4) ^ n)
        = Real.log ‖riemannZeta c‖ - n * Real.log (R / 4) := by
      rw [Real.log_div hζcpos.ne' (by positivity), Real.log_pow]
    have h4 : Real.log (3 * R / 4) = Real.log 3 + Real.log (R / 4) := by
      rw [show 3 * R / 4 = 3 * (R / 4) by ring, Real.log_mul (by norm_num) (by positivity)]
    have h5 : L = Real.log (M₀ + 1) - Real.log ‖riemannZeta c‖ := by
      rw [hLdef, Real.log_div hMpos.ne' hζcpos.ne']
    have h6 : (0 : ℝ) ≤ n * Real.log 3 := by
      have := Real.log_pos (by norm_num : (1 : ℝ) < 3)
      positivity
    have h4n : (n : ℝ) * Real.log (3 * R / 4)
        = (n : ℝ) * Real.log 3 + (n : ℝ) * Real.log (R / 4) := by
      rw [h4]; ring
    rw [h1, h5]
    have := h3 ▸ h2
    linarith
  -- assemble
  calc ‖logDeriv riemannZeta z‖
      ≤ (Real.log ((M₀ + 1) / ‖riemannZeta c‖) / Real.log 3) / δ
        + 32 * (Real.log ((M₀ + 1) / (3 * R / 4) ^ n) - Real.log ‖g c‖) / R := hmaster
    _ ≤ (L / Real.log 3) / δ + 32 * L / R := by
        gcongr
    _ = L * (1 / (δ * Real.log 3) + 32 / R) := by
        ring

#print axioms norm_logDeriv_riemannZeta_le_of_ball
