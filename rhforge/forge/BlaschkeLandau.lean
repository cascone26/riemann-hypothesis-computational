import BlaschkeDivide
import SharpBorelCaratheodory

/-!
# S3028 — one-radius Blaschke–Landau lemma (lane 1, κ 2 → 1 program, brick 4a: generic, centre 0)

For F with F(w) = (∏_{a∈S}(w − a)^{m a})·g₀(w) on ball 0 R' (R < R'), S ⊆ ball 0 R (the zeros inside the disk), g₀ analytic
on ball 0 R' and zero-free on ball 0 R, ‖F‖ ≤ Ms < Mζ on the sphere, F(0) ≠ 0, Re a ≤ 0 on S and a₀ ∈ S:

  Re(−F'/F(0)) ≤ 2·(log Mζ − log‖F 0‖)/R − m a₀·(Re(1/(−a₀)) + Re a₀/R²).

Proof: F = B·G (S3026) with G zero-free on the WHOLE ball 0 R, ‖G‖ = ‖F‖ ≤ Ms on the sphere ⇒ (max modulus) ‖G‖ ≤ Ms on the
disk; sharp BC (S3018) on radius R for G; Blaschke part by S3024; ‖G 0‖ ≥ ‖F 0‖ since ‖B 0‖ ≤ 1.
Compare the two-radii pipeline: 2·log(…)/(R/2) — this is HALF the growth term (κ 2 → 1). Zeros ON the sphere are harmless
(they stay in g₀; only ‖G‖ = ‖F‖ on the sphere is used).
-/

open Complex Metric Finset

namespace Blaschke

theorem landau_blaschke {F g₀ : ℂ → ℂ} {R R' Ms Mζ : ℝ} {S : Finset ℂ} (m : ℂ → ℕ) {a₀ : ℂ}
    (hR : 0 < R) (hRR' : R < R')
    (hS : ∀ a ∈ S, ‖a‖ < R) (h0 : (0 : ℂ) ∉ S) (hre : ∀ a ∈ S, a.re ≤ 0) (ha₀ : a₀ ∈ S)
    (hg₀ : AnalyticOnNhd ℂ g₀ (ball 0 R')) (hg₀ne : ∀ w ∈ ball (0 : ℂ) R, g₀ w ≠ 0)
    (hfac : ∀ w ∈ ball (0 : ℂ) R', F w = (∏ a ∈ S, (w - a) ^ m a) * g₀ w)
    (hM : ∀ w ∈ sphere (0 : ℂ) R, ‖F w‖ ≤ Ms) (hsl : Ms < Mζ) (hF0 : F 0 ≠ 0) :
    (-logDeriv F 0).re
      ≤ 2 * (Real.log Mζ - Real.log ‖F 0‖) / R - (m a₀ : ℝ) * ((1 / (-a₀)).re + a₀.re / R ^ 2) := by
  set G := Gfac R S m g₀ with hGdef
  set B := bprod R S m with hBdef
  have hcl : ∀ w ∈ closedBall (0 : ℂ) R, w ∈ ball (0 : ℂ) R' := fun w hw =>
    closedBall_subset_ball hRR' hw
  have hnormw : ∀ w ∈ closedBall (0 : ℂ) R, ‖w‖ ≤ R := fun w hw => by simpa using hw
  have hFBG : ∀ w ∈ closedBall (0 : ℂ) R, F w = B w * G w := fun w hw => by
    rw [hfac w (hcl w hw)]; exact fac_eq_bprod_mul m g₀ hS (hnormw w hw)
  -- G analytic on ball 0 R'
  have hGan : AnalyticOnNhd ℂ G (ball 0 R') := by
    have hpoly : AnalyticOnNhd ℂ (fun w : ℂ => ∏ a ∈ S, (((R : ℂ) ^ 2 - (starRingEnd ℂ) a * w) / (R : ℂ)) ^ m a)
        (ball 0 R') := by
      have heq : (fun w : ℂ => ∏ a ∈ S, (((R : ℂ) ^ 2 - (starRingEnd ℂ) a * w) / (R : ℂ)) ^ m a)
          = ∏ a ∈ S, (fun w : ℂ => (((R : ℂ) ^ 2 - (starRingEnd ℂ) a * w) / (R : ℂ)) ^ m a) := by
        ext w; simp [Finset.prod_apply]
      rw [heq]
      apply Finset.analyticOnNhd_prod
      intro a _ z _
      exact ((analyticAt_const.sub (analyticAt_const.mul analyticAt_id)).div_const).pow _
    exact hg₀.mul hpoly
  have hGdiff : DifferentiableOn ℂ G (ball 0 R) := (hGan.mono (ball_subset_ball hRR'.le)).differentiableOn
  have hGne : ∀ w ∈ ball (0 : ℂ) R, G w ≠ 0 := fun w hw =>
    Gfac_ne_zero m hS (hnormw w (ball_subset_closedBall hw)) (hg₀ne w hw)
  -- max modulus: ‖G‖ ≤ Ms on the disk
  have hGsphere : ∀ w ∈ sphere (0 : ℂ) R, ‖G w‖ ≤ Ms := by
    intro w hw
    have hwn : ‖w‖ = R := by simpa using hw
    have := Gfac_norm_sphere m g₀ hS hwn
    rw [← hfac w (hcl w (sphere_subset_closedBall hw))] at this
    rw [hGdef, this]; exact hM w hw
  have hGball : ∀ w ∈ ball (0 : ℂ) R, ‖G w‖ ≤ Ms := by
    intro w hw
    have hd : DiffContOnCl ℂ G (ball (0 : ℂ) R) :=
      ⟨hGdiff, by rw [closure_ball _ hR.ne']; exact (hGan.mono (closedBall_subset_ball hRR')).continuousOn⟩
    refine Complex.norm_le_of_forall_mem_frontier_norm_le isBounded_ball hd ?_ (subset_closure hw)
    rw [frontier_ball _ hR.ne']; exact hGsphere
  have h0ball : (0 : ℂ) ∈ ball (0 : ℂ) R := mem_ball_self hR
  have hG0 : G 0 ≠ 0 := hGne 0 h0ball
  have hB0 : B 0 ≠ 0 := bprod_ne_zero m hS (by simp [hR.le]) h0
  have hGlt : ‖G 0‖ < Mζ := lt_of_le_of_lt (hGball 0 h0ball) hsl
  have hbc := SharpBC.re_logDeriv_ge_of_bounded_sharp hR hGdiff hGne
    (fun w hw => (hGball w hw).trans hsl.le) hGlt
  -- logDeriv split
  have hev : F =ᶠ[nhds 0] fun w => B w * G w := by
    filter_upwards [isOpen_ball.mem_nhds h0ball] with w hw
    exact hFBG w (ball_subset_closedBall hw)
  have hBd : DifferentiableAt ℂ B 0 := by
    have : B = fun w => ∏ a ∈ S, b R a w ^ m a := rfl
    rw [this]
    exact DifferentiableAt.fun_finsetProd (fun a ha =>
      (b_differentiableAt (hS a ha) (by simp [(lt_of_le_of_lt (norm_nonneg _) (hS a ha)).le])).pow _)
  have hGd : DifferentiableAt ℂ G 0 := (hGan 0 (mem_ball_self (hR.trans hRR'))).differentiableAt
  have hsplit := logDeriv_split hev hB0 hG0 hBd hGd
  have hBpart := bprod_neg_logDeriv_re_le m hS h0 hre ha₀
  -- ‖F 0‖ ≤ ‖G 0‖
  have hF0G : ‖F 0‖ ≤ ‖G 0‖ := by
    rw [hFBG 0 (mem_closedBall_self hR.le), norm_mul]
    have := bprod_norm_le_one m hS (w := 0) (by simp [hR.le])
    nlinarith [norm_nonneg (G 0), norm_nonneg (B 0)]
  have hlog : Real.log ‖F 0‖ ≤ Real.log ‖G 0‖ := Real.log_le_log (norm_pos_iff.mpr hF0) hF0G
  have hre_eq : (-logDeriv F 0).re = (-logDeriv B 0).re + (-logDeriv G 0).re := by
    rw [hsplit, neg_add, Complex.add_re]
  rw [hre_eq, Complex.neg_re (logDeriv G 0)]
  have h2R : 2 * (Real.log Mζ - Real.log ‖G 0‖) / R ≤ 2 * (Real.log Mζ - Real.log ‖F 0‖) / R := by
    apply div_le_div_of_nonneg_right _ hR.le; linarith
  linarith

end Blaschke

#print axioms Blaschke.landau_blaschke
