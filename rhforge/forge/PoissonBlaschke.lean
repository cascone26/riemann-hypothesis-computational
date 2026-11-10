import BlaschkeLandau
import SchwarzDeriv
import ExistsHoloBranchBall
import ReEqLogNormOfExp

/-!
# S3039 — Poisson–Blaschke identity/bound (lane 1, Poisson–Landau brick P2b, generic centre 0)

F = (∏_{a∈S}(w−a)^{m a})·g₀ on ball 0 R', S ⊆ ball 0 r, g₀ analytic on ball 0 R' and zero-free on ball 0 r' (r < r' ≤ R'),
0 ∉ S, Re a ≤ 0 on S, a₀ ∈ S, F(0) ≠ 0:

  Re(−F'/F(0)) ≤ −avg_{|ζ|=r} Re((2ζ/ζ²)·log‖F ζ‖) − m a₀·(Re(1/(−a₀)) + Re a₀/r²).

i.e. the growth term is the CIRCLE AVERAGE of (2Re ζ/r²)·log|F| — not its maximum (contrast S3028).
Proof: G = Gfac r S m g₀ (S3026) is zero-free on ball 0 ρ, ρ = min(r', R', r²/μ') > r (μ' ≥ max‖a‖, μ' < r);
log branch f of G on ball 0 ρ; Schwarz formula (S3034) at radius r for f; Re f = log‖G‖ = log‖F‖ on the sphere; logDeriv
split F = B·G; Blaschke part (S3024).
-/

open Complex Metric Finset

namespace PoissonBlaschke

theorem exists_norm_bound (S : Finset ℂ) {r : ℝ} (hr : 0 < r) (hS : ∀ a ∈ S, ‖a‖ < r) :
    ∃ μ : ℝ, 0 < μ ∧ μ < r ∧ ∀ a ∈ S, ‖a‖ ≤ μ := by
  induction S using Finset.induction_on with
  | empty => exact ⟨r / 2, by linarith, by linarith, by simp⟩
  | insert b T hb ih =>
    obtain ⟨μ, h1, h2, h3⟩ := ih (fun a ha => hS a (Finset.mem_insert_of_mem ha))
    refine ⟨max μ ‖b‖, lt_max_of_lt_left h1, max_lt h2 (hS b (Finset.mem_insert_self b T)), ?_⟩
    intro a ha
    rcases Finset.mem_insert.mp ha with rfl | ha
    · exact le_max_right _ _
    · exact (h3 a ha).trans (le_max_left _ _)

theorem poisson_blaschke {F g₀ : ℂ → ℂ} {r r' R' : ℝ} {S : Finset ℂ} (m : ℂ → ℕ) {a₀ : ℂ}
    (hr : 0 < r) (hrr' : r < r') (hr'R' : r' ≤ R')
    (hS : ∀ a ∈ S, ‖a‖ < r) (h0 : (0 : ℂ) ∉ S) (hre : ∀ a ∈ S, a.re ≤ 0) (ha₀ : a₀ ∈ S)
    (hg₀ : AnalyticOnNhd ℂ g₀ (ball 0 R')) (hg₀ne : ∀ w ∈ ball (0 : ℂ) r', g₀ w ≠ 0)
    (hfac : ∀ w ∈ ball (0 : ℂ) R', F w = (∏ a ∈ S, (w - a) ^ m a) * g₀ w) (hF0 : F 0 ≠ 0) :
    (-logDeriv F 0).re
      ≤ -Real.circleAverage (fun ζ => ((2 * ζ / (ζ - 0) ^ 2) * ((Real.log ‖F ζ‖ : ℝ) : ℂ)).re) 0 r
        - (m a₀ : ℝ) * ((1 / (-a₀)).re + a₀.re / r ^ 2) := by
  set G := Blaschke.Gfac r S m g₀ with hGdef
  set B := Blaschke.bprod r S m with hBdef
  obtain ⟨μ, hμ0, hμr, hμ⟩ := exists_norm_bound S hr hS
  set ρ := min r' (r ^ 2 / μ) with hρ
  have hrμ : r < r ^ 2 / μ := by rw [lt_div_iff₀ hμ0]; nlinarith
  have hrρ : r < ρ := lt_min hrr' hrμ
  have hρpos : 0 < ρ := hr.trans hrρ
  have hρR' : ρ ≤ R' := (min_le_left _ _).trans hr'R'
  -- G analytic on ball 0 R'
  have hGan : AnalyticOnNhd ℂ G (ball 0 R') := by
    have hpoly : AnalyticOnNhd ℂ (fun w : ℂ => ∏ a ∈ S, (((r : ℂ) ^ 2 - (starRingEnd ℂ) a * w) / (r : ℂ)) ^ m a)
        (ball 0 R') := by
      have heq : (fun w : ℂ => ∏ a ∈ S, (((r : ℂ) ^ 2 - (starRingEnd ℂ) a * w) / (r : ℂ)) ^ m a)
          = ∏ a ∈ S, (fun w : ℂ => (((r : ℂ) ^ 2 - (starRingEnd ℂ) a * w) / (r : ℂ)) ^ m a) := by
        ext w; simp [Finset.prod_apply]
      rw [heq]
      apply Finset.analyticOnNhd_prod
      intro a _ z _
      exact ((analyticAt_const.sub (analyticAt_const.mul analyticAt_id)).div_const).pow _
    exact hg₀.mul hpoly
  -- G zero-free on ball 0 ρ
  have hGne : ∀ w ∈ ball (0 : ℂ) ρ, G w ≠ 0 := by
    intro w hw
    have hwn : ‖w‖ < ρ := by simpa using hw
    have hwr' : w ∈ ball (0 : ℂ) r' := by
      rw [mem_ball, dist_zero_right]; exact lt_of_lt_of_le hwn (min_le_left _ _)
    have hwμ : ‖w‖ < r ^ 2 / μ := lt_of_lt_of_le hwn (min_le_right _ _)
    simp only [hGdef, Blaschke.Gfac]
    refine mul_ne_zero (hg₀ne w hwr') (Finset.prod_ne_zero_iff.mpr (fun a ha => pow_ne_zero _ ?_))
    refine div_ne_zero ?_ (by exact_mod_cast hr.ne')
    intro h
    have h1 : ‖(starRingEnd ℂ) a * w‖ = ‖((r : ℂ) ^ 2)‖ := by rw [sub_eq_zero.mp h]
    rw [norm_mul, Complex.norm_conj, norm_pow, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hr] at h1
    have h2 : ‖a‖ * ‖w‖ < r ^ 2 := by
      calc ‖a‖ * ‖w‖ ≤ μ * ‖w‖ := mul_le_mul_of_nonneg_right (hμ a ha) (norm_nonneg _)
        _ < μ * (r ^ 2 / μ) := mul_lt_mul_of_pos_left hwμ hμ0
        _ = r ^ 2 := by field_simp
    linarith
  have hGdiffρ : DifferentiableOn ℂ G (ball 0 ρ) :=
    (hGan.mono (ball_subset_ball hρR')).differentiableOn
  -- log branch on ball 0 ρ and Schwarz at radius r
  obtain ⟨f, hfd, hfe⟩ := exists_holo_branch_on_ball hρpos hGdiffρ hGne
  have hfdiff : DifferentiableOn ℂ f (ball 0 ρ) := fun z hz => (hfd z hz).differentiableAt.differentiableWithinAt
  have hfan : AnalyticOnNhd ℂ f (ball 0 ρ) := hfdiff.analyticOnNhd isOpen_ball
  have hsch := SchwarzDeriv.deriv_eq_schwarz hr hrρ hfan
  have h0ρ : (0 : ℂ) ∈ ball (0 : ℂ) ρ := mem_ball_self hρpos
  have hdf0 : deriv f 0 = logDeriv G 0 := (hfd 0 h0ρ).deriv
  -- real part of the circle average
  have hcont : ContinuousOn (fun ζ : ℂ => (2 * ζ / (ζ - 0) ^ 2) • (((f ζ).re : ℝ) : ℂ)) (sphere (0 : ℂ) |r|) := by
    rw [abs_of_pos hr]
    have hfc : ContinuousOn f (sphere (0 : ℂ) r) :=
      hfan.continuousOn.mono (sphere_subset_closedBall.trans (closedBall_subset_ball hrρ))
    show ContinuousOn (fun ζ : ℂ => (2 * ζ / (ζ - 0) ^ 2) * (((f ζ).re : ℝ) : ℂ)) _
    refine ContinuousOn.mul ?_ (Complex.continuous_ofReal.comp_continuousOn
      (Complex.continuous_re.comp_continuousOn hfc))
    refine ContinuousOn.div (by fun_prop) (by fun_prop) ?_
    intro ζ hζ
    have : ζ ≠ 0 := by
      intro h; rw [h, mem_sphere, dist_self] at hζ; linarith
    simpa using this
  have hci := hcont.circleIntegrable'
  have hre_avg := Complex.reCLM.restrictScalars ℝ |>.circleAverage_comp_comm hci
  have hGre : (logDeriv G 0).re
      = Real.circleAverage (fun ζ => ((2 * ζ / (ζ - 0) ^ 2) * ((Real.log ‖F ζ‖ : ℝ) : ℂ)).re) 0 r := by
    rw [← hdf0, hsch]
    have : (Real.circleAverage (fun ζ => (2 * ζ / (ζ - 0) ^ 2) • (((f ζ).re : ℝ) : ℂ)) 0 r).re
        = Real.circleAverage (fun ζ => ((2 * ζ / (ζ - 0) ^ 2) • (((f ζ).re : ℝ) : ℂ)).re) 0 r := by
      exact hre_avg.symm
    rw [this]
    apply Real.circleAverage_congr_sphere
    intro ζ hζ
    rw [abs_of_pos hr] at hζ
    have hζn : ‖ζ‖ = r := by simpa using hζ
    have hζρ : ζ ∈ ball (0 : ℂ) ρ := sphere_subset_closedBall.trans (closedBall_subset_ball hrρ) hζ
    have hlog : (f ζ).re = Real.log ‖F ζ‖ := by
      rw [re_eq_log_norm_of_exp (hfe ζ hζρ)]
      have := Blaschke.Gfac_norm_sphere m g₀ hS hζn
      rw [← hfac ζ (sphere_subset_closedBall.trans (closedBall_subset_ball (hrr'.trans_le hr'R')) hζ)] at this
      rw [hGdef, this]
    simp only [smul_eq_mul, hlog]
  -- logDeriv split F = B·G near 0
  have hnormw : ∀ w ∈ closedBall (0 : ℂ) r, ‖w‖ ≤ r := fun w hw => by simpa using hw
  have hFBG : ∀ w ∈ closedBall (0 : ℂ) r, F w = B w * G w := fun w hw => by
    rw [hfac w (closedBall_subset_ball (hrr'.trans_le hr'R') hw)]
    exact Blaschke.fac_eq_bprod_mul m g₀ hS (hnormw w hw)
  have h0r : (0 : ℂ) ∈ ball (0 : ℂ) r := mem_ball_self hr
  have hev : F =ᶠ[nhds 0] fun w => B w * G w := by
    filter_upwards [isOpen_ball.mem_nhds h0r] with w hw
    exact hFBG w (ball_subset_closedBall hw)
  have hB0 : B 0 ≠ 0 := Blaschke.bprod_ne_zero m hS (by simp [hr.le]) h0
  have hG0 : G 0 ≠ 0 := hGne 0 h0ρ
  have hBd : DifferentiableAt ℂ B 0 := by
    have : B = fun w => ∏ a ∈ S, Blaschke.b r a w ^ m a := rfl
    rw [this]
    exact DifferentiableAt.fun_finsetProd (fun a ha =>
      (Blaschke.b_differentiableAt (hS a ha) (by simp [(lt_of_le_of_lt (norm_nonneg _) (hS a ha)).le])).pow _)
  have hGd : DifferentiableAt ℂ G 0 := (hGan 0 (mem_ball_self (hr.trans (hrr'.trans_le hr'R')))).differentiableAt
  have hsplit := Blaschke.logDeriv_split hev hB0 hG0 hBd hGd
  have hBpart := Blaschke.bprod_neg_logDeriv_re_le m hS h0 hre ha₀
  have hre_eq : (-logDeriv F 0).re = (-logDeriv B 0).re + (-logDeriv G 0).re := by
    rw [hsplit, neg_add, Complex.add_re]
  rw [hre_eq, Complex.neg_re (logDeriv G 0), hGre]
  linarith

end PoissonBlaschke

#print axioms PoissonBlaschke.exists_norm_bound
#print axioms PoissonBlaschke.poisson_blaschke
