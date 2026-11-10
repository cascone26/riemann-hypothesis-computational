import PoissonBlaschke

/-!
# S3042 — Poisson–Blaschke bound with ALL zeros dropped (lane 1, Poisson–Landau P5a)

`bprod_neg_logDeriv_re_nonpos`: Re a ≤ 0 on S ⇒ Re(−B'/B(0)) ≤ 0 (every Blaschke bracket is ≥ 0, S3022).
`poisson_blaschke_drop`: same hypotheses as S3039 without a tracked zero (S may be empty) ⇒
   Re(−F'/F(0)) ≤ −avg_{|ζ|=r} Re((2ζ/ζ²)·log‖F ζ‖)
— the growth form used for the k = 2, 3, 4 terms of the quartic region.
-/

open Complex Metric Finset

namespace PoissonBlaschke

theorem bprod_neg_logDeriv_re_nonpos {r : ℝ} {S : Finset ℂ} (m : ℂ → ℕ) (hS : ∀ a ∈ S, ‖a‖ < r)
    (h0 : (0 : ℂ) ∉ S) (hre : ∀ a ∈ S, a.re ≤ 0) (hr : 0 < r) :
    (-logDeriv (Blaschke.bprod r S m) 0).re ≤ 0 := by
  rw [Blaschke.bprod_logDeriv_zero m hS h0, Complex.neg_re, Complex.re_sum]
  have : ∀ a ∈ S, 0 ≤ ((m a : ℂ) * (-1 / a + (starRingEnd ℂ) a / (r : ℂ) ^ 2)).re := by
    intro a ha
    have e1 : (-1 / a : ℂ) = 1 / (-a) := by rw [div_neg, neg_div]
    have e2 : (r : ℂ) ^ 2 = ((r ^ 2 : ℝ) : ℂ) := by push_cast; ring
    have e3 : (m a : ℂ) = ((m a : ℝ) : ℂ) := by push_cast; ring
    rw [e1, e2, e3, Complex.re_ofReal_mul, Complex.add_re, Complex.div_ofReal_re, Complex.conj_re]
    exact mul_nonneg (Nat.cast_nonneg _)
      (Blaschke.blaschke_zero_term_nonneg hr (hS a ha).le (fun h => h0 (h ▸ ha)) (hre a ha))
  have := Finset.sum_nonneg this
  linarith

theorem poisson_blaschke_drop {F g₀ : ℂ → ℂ} {r r' R' : ℝ} {S : Finset ℂ} (m : ℂ → ℕ)
    (hr : 0 < r) (hrr' : r < r') (hr'R' : r' ≤ R')
    (hS : ∀ a ∈ S, ‖a‖ < r) (h0 : (0 : ℂ) ∉ S) (hre : ∀ a ∈ S, a.re ≤ 0)
    (hg₀ : AnalyticOnNhd ℂ g₀ (ball 0 R')) (hg₀ne : ∀ w ∈ ball (0 : ℂ) r', g₀ w ≠ 0)
    (hfac : ∀ w ∈ ball (0 : ℂ) R', F w = (∏ a ∈ S, (w - a) ^ m a) * g₀ w) (hF0 : F 0 ≠ 0) :
    (-logDeriv F 0).re
      ≤ -Real.circleAverage (fun ζ => ((2 * ζ / (ζ - 0) ^ 2) * ((Real.log ‖F ζ‖ : ℝ) : ℂ)).re) 0 r := by
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
  have hBpart := bprod_neg_logDeriv_re_nonpos m hS h0 hre hr
  have hre_eq : (-logDeriv F 0).re = (-logDeriv B 0).re + (-logDeriv G 0).re := by
    rw [hsplit, neg_add, Complex.add_re]
  rw [hre_eq, Complex.neg_re (logDeriv G 0), hGre]
  linarith

end PoissonBlaschke

#print axioms PoissonBlaschke.bprod_neg_logDeriv_re_nonpos
#print axioms PoissonBlaschke.poisson_blaschke_drop
