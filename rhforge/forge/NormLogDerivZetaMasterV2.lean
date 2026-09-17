import NormLogDerivQuarterBall
import ZeroCountLog
import LogDerivProdPow

set_option maxHeartbeats 1000000

open Complex Metric

/-!
# §1680 ★★★ — NORM-LANDAU MASTER v2 (consistent quarter-ball geometry, the ζ-usable form):
`‖ζ'/ζ(z)‖ ≤ (log(M/‖ζ c‖)/log 3)/δ + 32(log(M/(3R/4)^Σm) − log‖g c‖)/R`

§1677 rewired to §1679's ζ-compatible parameterization: `S` = zeros in the quarter-ball (matching
the §1674 count REGION exactly), `g` nonvanishing on the quarter-ball only (satisfiable for ζ),
bound valid for `z ∈ closedBall c (R/16)` keeping `δ` from every zero of `S`. This is the form the
box-edge instantiation consumes: `c = 1 + 1/log T + iT`, `R = 1/4`, `M = O(T)`,
`‖ζ c‖ ≥ 1/O(log T)` (§1678), `δ = c₀/(2 log(T+2))` (§1667's zero-free geometry) ⟹
`‖ζ'/ζ‖ = O(log²T)` on the connectors.
-/

theorem norm_logDeriv_riemannZeta_le_master' {c z : ℂ} {R M δ : ℝ} {g P : ℂ → ℂ}
    {S : Finset ℂ} {m : ℂ → ℕ}
    (hR : 0 < R) (hδ : 0 < δ)
    (hPdef : P = fun s => ∏ ρ ∈ S, (s - ρ) ^ (m ρ))
    (hS : ∀ ρ ∈ S, ρ ∈ closedBall c (R / 4))
    (hfac : ∀ s ∈ closedBall c R, riemannZeta s = P s * g s)
    (hζsphere : ∀ s ∈ sphere c R, ‖riemannZeta s‖ ≤ M)
    (hgdiff : DiffContOnCl ℂ g (ball c R))
    (hgnz : ∀ u ∈ ball c (R / 4), g u ≠ 0)
    (hgc_lt : ‖g c‖ < M / (3 * R / 4) ^ (∑ ρ ∈ S, m ρ))
    (hζc0 : riemannZeta c ≠ 0)
    (hz : z ∈ closedBall c (R / 16))
    (hdist : ∀ ρ ∈ S, δ ≤ ‖z - ρ‖) :
    ‖logDeriv riemannZeta z‖
      ≤ (Real.log (M / ‖riemannZeta c‖) / Real.log 3) / δ
        + 32 * (Real.log (M / (3 * R / 4) ^ (∑ ρ ∈ S, m ρ)) - Real.log ‖g c‖) / R := by
  have hzcb : z ∈ closedBall c R := by
    rw [mem_closedBall] at hz ⊢
    linarith
  have hzball : z ∈ ball c R := by
    rw [mem_closedBall] at hz
    rw [mem_ball]
    linarith
  have hzquarter : z ∈ ball c (R / 4) := by
    rw [mem_closedBall] at hz
    rw [mem_ball]
    linarith
  have hsubne : ∀ ρ ∈ S, z - ρ ≠ 0 := by
    intro ρ hρ h0
    have := hdist ρ hρ
    rw [h0, norm_zero] at this
    linarith
  have hPz : P z ≠ 0 := by
    rw [hPdef]
    exact Finset.prod_ne_zero_iff.mpr (fun ρ hρ => pow_ne_zero _ (hsubne ρ hρ))
  have hgz : g z ≠ 0 := hgnz z hzquarter
  have hPd : DifferentiableAt ℂ P z := by
    rw [hPdef]
    apply DifferentiableAt.fun_finsetProd
    intro ρ _
    exact (differentiableAt_id.sub_const ρ).pow _
  have hgd : DifferentiableAt ℂ g z :=
    hgdiff.differentiableOn.differentiableAt (isOpen_ball.mem_nhds hzball)
  have hev : riemannZeta =ᶠ[nhds z] fun s => P s * g s := by
    filter_upwards [isOpen_ball.mem_nhds hzball] with s hs
    exact hfac s (ball_subset_closedBall hs)
  have hsplit : logDeriv riemannZeta z = logDeriv P z + logDeriv g z := by
    have h1 : logDeriv riemannZeta z = logDeriv (fun s => P s * g s) z := by
      rw [logDeriv_apply, logDeriv_apply, hev.deriv_eq, hfac z hzcb]
    rw [h1, logDeriv_mul z hPz hgz hPd hgd]
  have hcount := sum_mult_le_log_of_factor hR hPdef hS
    (fun s hs => hfac s (sphere_subset_closedBall hs))
    (hfac c (mem_closedBall_self (by linarith)))
    hζsphere hgdiff hζc0
  have hPterm : ‖logDeriv P z‖ ≤ (Real.log (M / ‖riemannZeta c‖) / Real.log 3) / δ := by
    rw [hPdef, logDeriv_prod_pow hsubne]
    calc ‖∑ ρ ∈ S, (m ρ : ℂ) / (z - ρ)‖
        ≤ ∑ ρ ∈ S, ‖(m ρ : ℂ) / (z - ρ)‖ := norm_sum_le _ _
      _ ≤ ∑ ρ ∈ S, (m ρ : ℝ) / δ := by
          apply Finset.sum_le_sum
          intro ρ hρ
          rw [norm_div]
          have hnum : ‖((m ρ : ℕ) : ℂ)‖ = (m ρ : ℝ) := by
            rw [show ((m ρ : ℕ) : ℂ) = ((m ρ : ℝ) : ℂ) by push_cast; ring,
              Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (Nat.cast_nonneg _)]
          rw [hnum]
          exact div_le_div_of_nonneg_left (Nat.cast_nonneg _) hδ (hdist ρ hρ)
      _ = ((∑ ρ ∈ S, m ρ : ℕ) : ℝ) / δ := by
          rw [← Finset.sum_div]
          congr 1
          push_cast
          ring
      _ ≤ (Real.log (M / ‖riemannZeta c‖) / Real.log 3) / δ := by
          gcongr
  have hgterm := norm_logDeriv_g_le_of_quarterBall_factor hR hPdef hS
    (fun s hs => hfac s (sphere_subset_closedBall hs)) hζsphere hgdiff hgnz hgc_lt z hz
  calc ‖logDeriv riemannZeta z‖ = ‖logDeriv P z + logDeriv g z‖ := by rw [hsplit]
    _ ≤ ‖logDeriv P z‖ + ‖logDeriv g z‖ := norm_add_le _ _
    _ ≤ (Real.log (M / ‖riemannZeta c‖) / Real.log 3) / δ
        + 32 * (Real.log (M / (3 * R / 4) ^ (∑ ρ ∈ S, m ρ)) - Real.log ‖g c‖) / R := by
          gcongr

#print axioms norm_logDeriv_riemannZeta_le_master'
