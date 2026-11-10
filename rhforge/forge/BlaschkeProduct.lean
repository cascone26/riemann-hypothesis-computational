import BlaschkeFactor

/-!
# S3024 — finite Blaschke product (lane 1, κ 2 → 1 program, brick 2)

B(w) = ∏_{a ∈ S} b_R(a, w)^{m a} over a finite set S of (centred) zeros with |a| < R:
* `bprod_norm_sphere`  — ‖w‖ = R ⇒ ‖B w‖ = 1   (so |ζ/B| = |ζ| on the sphere);
* `bprod_norm_le_one`  — ‖w‖ ≤ R ⇒ ‖B w‖ ≤ 1;
* `bprod_ne_zero`      — ‖w‖ ≤ R, w ∉ S ⇒ B w ≠ 0;
* `bprod_logDeriv_zero` — 0 ∉ S ⇒ logDeriv B 0 = Σ m a · (−1/a + ā/R²);
* `bprod_neg_logDeriv_re_le` — all a ∈ S with Re a ≤ 0: Re(−logDeriv B 0) ≤ −m a₀·(Re(1/(−a₀)) + Re a₀/R²) for any a₀ ∈ S
  (the Blaschke-corrected zero sum: drop every zero but a₀, each dropped bracket is ≥ 0 by S3022).
-/

open Complex Finset

namespace Blaschke

noncomputable def bprod (R : ℝ) (S : Finset ℂ) (m : ℂ → ℕ) (w : ℂ) : ℂ := ∏ a ∈ S, b R a w ^ m a

theorem b_ne_zero {R : ℝ} {a w : ℂ} (ha : ‖a‖ < R) (hw : ‖w‖ ≤ R) (hwa : w ≠ a) : b R a w ≠ 0 := by
  have hR : 0 < R := lt_of_le_of_lt (norm_nonneg _) ha
  unfold b
  exact div_ne_zero (mul_ne_zero (by exact_mod_cast hR.ne') (sub_ne_zero.mpr hwa)) (blaschke_den_ne ha hw)

theorem b_differentiableAt {R : ℝ} {a w : ℂ} (ha : ‖a‖ < R) (hw : ‖w‖ ≤ R) : DifferentiableAt ℂ (b R a) w := by
  have hden := blaschke_den_ne ha hw
  have : b R a = fun w => (R : ℂ) * (w - a) / ((R : ℂ) ^ 2 - (starRingEnd ℂ) a * w) := rfl
  rw [this]
  exact ((differentiableAt_id.sub_const a).const_mul _).div
    ((differentiableAt_id.const_mul _).const_sub _) hden

theorem bprod_norm_sphere {R : ℝ} {S : Finset ℂ} (m : ℂ → ℕ) (hS : ∀ a ∈ S, ‖a‖ < R) {w : ℂ} (hw : ‖w‖ = R) :
    ‖bprod R S m w‖ = 1 := by
  unfold bprod
  rw [norm_prod]
  exact Finset.prod_eq_one (fun a ha => by rw [norm_pow, blaschke_norm_sphere (hS a ha) hw, one_pow])

theorem bprod_norm_le_one {R : ℝ} {S : Finset ℂ} (m : ℂ → ℕ) (hS : ∀ a ∈ S, ‖a‖ < R) {w : ℂ} (hw : ‖w‖ ≤ R) :
    ‖bprod R S m w‖ ≤ 1 := by
  unfold bprod
  rw [norm_prod]
  exact Finset.prod_le_one (fun a _ => norm_nonneg _)
    (fun a ha => by rw [norm_pow]; exact pow_le_one₀ (norm_nonneg _) (blaschke_norm_le_one (hS a ha) hw))

theorem bprod_ne_zero {R : ℝ} {S : Finset ℂ} (m : ℂ → ℕ) (hS : ∀ a ∈ S, ‖a‖ < R) {w : ℂ} (hw : ‖w‖ ≤ R)
    (hwS : w ∉ S) : bprod R S m w ≠ 0 := by
  unfold bprod
  exact Finset.prod_ne_zero_iff.mpr (fun a ha => pow_ne_zero _ (b_ne_zero (hS a ha) hw (fun h => hwS (h ▸ ha))))

theorem bprod_logDeriv_zero {R : ℝ} {S : Finset ℂ} (m : ℂ → ℕ) (hS : ∀ a ∈ S, ‖a‖ < R) (h0 : (0 : ℂ) ∉ S) :
    logDeriv (bprod R S m) 0 = ∑ a ∈ S, (m a : ℂ) * (-1 / a + (starRingEnd ℂ) a / (R : ℂ) ^ 2) := by
  have hR0 : ∀ a ∈ S, ‖(0 : ℂ)‖ ≤ R := fun a ha => by
    rw [norm_zero]; exact le_trans (norm_nonneg _) (hS a ha).le
  have hfun : bprod R S m = fun w => ∏ a ∈ S, (fun w => b R a w ^ m a) w := rfl
  rw [hfun, logDeriv_prod]
  · refine Finset.sum_congr rfl (fun a ha => ?_)
    rw [logDeriv_fun_pow (b_differentiableAt (hS a ha) (hR0 a ha)),
      blaschke_logDeriv_zero (hS a ha) (fun h => h0 (h ▸ ha))]
  · intro a ha
    exact pow_ne_zero _ (b_ne_zero (hS a ha) (hR0 a ha) (fun h => h0 (h ▸ ha)))
  · intro a ha
    exact (b_differentiableAt (hS a ha) (hR0 a ha)).pow _

theorem bprod_neg_logDeriv_re_le {R : ℝ} {S : Finset ℂ} (m : ℂ → ℕ) (hS : ∀ a ∈ S, ‖a‖ < R) (h0 : (0 : ℂ) ∉ S)
    (hre : ∀ a ∈ S, a.re ≤ 0) {a₀ : ℂ} (ha₀ : a₀ ∈ S) :
    (-logDeriv (bprod R S m) 0).re ≤ -((m a₀ : ℝ) * ((1 / (-a₀)).re + a₀.re / R ^ 2)) := by
  have hR : 0 < R := lt_of_le_of_lt (norm_nonneg _) (hS a₀ ha₀)
  rw [bprod_logDeriv_zero m hS h0]
  have hX : ∀ a ∈ S, 0 ≤ (m a : ℝ) * ((1 / (-a)).re + a.re / R ^ 2) := fun a ha =>
    mul_nonneg (Nat.cast_nonneg _)
      (blaschke_zero_term_nonneg hR (hS a ha).le (fun h => h0 (h ▸ ha)) (hre a ha))
  have hterm : ∀ a ∈ S, ((m a : ℂ) * (-1 / a + (starRingEnd ℂ) a / (R : ℂ) ^ 2)).re
      = (m a : ℝ) * ((1 / (-a)).re + a.re / R ^ 2) := by
    intro a _
    have e1 : (-1 / a : ℂ) = 1 / (-a) := by rw [div_neg, neg_div]
    have e2 : (R : ℂ) ^ 2 = ((R ^ 2 : ℝ) : ℂ) := by push_cast; ring
    have e3 : (m a : ℂ) = ((m a : ℝ) : ℂ) := by push_cast; ring
    rw [e1, e2, e3, Complex.re_ofReal_mul, Complex.add_re, Complex.div_ofReal_re, Complex.conj_re]
  rw [Complex.neg_re, Complex.re_sum, Finset.sum_congr rfl hterm]
  have := Finset.single_le_sum hX ha₀
  linarith

end Blaschke

#print axioms Blaschke.bprod_norm_sphere
#print axioms Blaschke.bprod_norm_le_one
#print axioms Blaschke.bprod_ne_zero
#print axioms Blaschke.bprod_logDeriv_zero
#print axioms Blaschke.bprod_neg_logDeriv_re_le
