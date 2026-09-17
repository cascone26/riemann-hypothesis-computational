import NormLogDerivBallCenter
import HgBFromFactorTwoRadii
import Mathlib.Analysis.Complex.AbsMax

set_option maxHeartbeats 1000000

open Complex Metric

/-!
# §1672 ★★ — NORM bound for the zero-free factor of `ζ`: `‖logDeriv g s‖ ≤ 16·log(M'/‖g c‖)/R`

The norm analogue of the tower's `re_logDeriv_ge_of_halfBall_factor'` (§HgBFromFactorTwoRadii),
upgrading its `Re ≥ -B` conclusion to a full `‖·‖ ≤` bound on the `R/8`-ball. Same hypotheses:
`ζ = P·g` on the sphere with `P` the zeros-in-the-half-ball polynomial, `‖ζ‖ ≤ Mζ` on the sphere.
Chain: the tower's sphere bound `‖g‖ ≤ M' = Mζ/(R/2)^Σm` (via `norm_prod_sub_half_ge` +
`norm_factor_le`) → MAX MODULUS (`norm_le_of_forall_mem_frontier_norm_le`) extends it to the ball
→ §1670's norm-Landau core on `ball c (R/2)` (where `g` is nonvanishing) with
`M'' = log(M'/‖g c‖)`. The last remaining Landau piece before the box-edge `O(log²T)`.
-/

theorem norm_logDeriv_g_le_of_halfBall_factor {c : ℂ} {R Mζ : ℝ} {g P : ℂ → ℂ} {S : Finset ℂ}
    {m : ℂ → ℕ} (hR : 0 < R) (hPdef : P = fun s => ∏ ρ ∈ S, (s - ρ) ^ (m ρ))
    (hS : ∀ ρ ∈ S, ρ ∈ closedBall c (R / 2))
    (hfac : ∀ s ∈ sphere c R, riemannZeta s = P s * g s)
    (hζsphere : ∀ s ∈ sphere c R, ‖riemannZeta s‖ ≤ Mζ)
    (hgdiff : DiffContOnCl ℂ g (ball c R)) (hgnz : ∀ z ∈ ball c (R / 2), g z ≠ 0)
    (hgc_lt : ‖g c‖ < Mζ / (R / 2) ^ (∑ ρ ∈ S, m ρ)) :
    ∀ s ∈ closedBall c (R / 8),
      ‖logDeriv g s‖ ≤ 16 * (Real.log (Mζ / (R / 2) ^ (∑ ρ ∈ S, m ρ)) - Real.log ‖g c‖) / R := by
  set M' : ℝ := Mζ / (R / 2) ^ (∑ ρ ∈ S, m ρ) with hM'def
  have hR2 : 0 < R / 2 := by linarith
  have hgc0 : g c ≠ 0 := hgnz c (mem_ball_self hR2)
  have hgcpos : 0 < ‖g c‖ := norm_pos_iff.mpr hgc0
  have hM'pos : 0 < M' := lt_trans hgcpos hgc_lt
  -- the tower's sphere bound for g
  have hgsphere : ∀ s ∈ sphere c R, ‖g s‖ ≤ M' := by
    intro s hs
    have hPlow : (R / 2) ^ (∑ ρ ∈ S, m ρ) ≤ ‖P s‖ := by
      rw [hPdef]
      exact norm_prod_sub_half_ge hR.le m hs hS
    exact norm_factor_le (hfac s hs) (by positivity) hPlow (hζsphere s hs)
  -- max modulus: extend to the closed ball
  have hgball : ∀ z ∈ ball c (R / 2), ‖g z‖ ≤ M' := by
    intro z hz
    have hfr : ∀ u ∈ frontier (ball c R), ‖g u‖ ≤ M' := by
      intro u hu
      rw [frontier_ball c hR.ne'] at hu
      exact hgsphere u hu
    have hcl : z ∈ closure (ball c R) := subset_closure (by
      rw [mem_ball] at hz ⊢
      linarith)
    exact Complex.norm_le_of_forall_mem_frontier_norm_le isBounded_ball hgdiff hfr hcl
  -- norm-Landau core on the half-ball
  set M'' : ℝ := Real.log (M' / ‖g c‖) with hM''def
  have hM''pos : 0 < M'' := Real.log_pos ((one_lt_div hgcpos).mpr hgc_lt)
  have hgrow : ∀ z ∈ ball c (R / 2), ‖g z‖ ≤ Real.exp M'' * ‖g c‖ := by
    intro z hz
    rw [hM''def, Real.exp_log (by positivity), div_mul_cancel₀ _ hgcpos.ne']
    exact hgball z hz
  have hgdiff2 : DifferentiableOn ℂ g (ball c (R / 2)) := by
    apply hgdiff.differentiableOn.mono
    intro z hz
    rw [mem_ball] at hz ⊢
    linarith
  have hcore := norm_logDeriv_le_of_forall_norm_le_center hR2 hM''pos hgdiff2 hgnz hgrow
  intro s hs
  have hs' : s ∈ closedBall c (R / 2 / 4) := by
    rw [mem_closedBall] at hs ⊢
    linarith
  calc ‖logDeriv g s‖ ≤ 8 * M'' / (R / 2) := hcore s hs'
    _ = 16 * M'' / R := by field_simp; ring
    _ = 16 * (Real.log M' - Real.log ‖g c‖) / R := by
        rw [hM''def, Real.log_div hM'pos.ne' hgcpos.ne']

#print axioms norm_logDeriv_g_le_of_halfBall_factor
