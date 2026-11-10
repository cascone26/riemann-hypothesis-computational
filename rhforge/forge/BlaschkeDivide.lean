import BlaschkeProduct

/-!
# S3026 — dividing by the Blaschke product (lane 1, κ 2 → 1 program, brick 3, generic form)

For f = (∏_{a∈S}(w − a)^{m a})·g₀ on the closed disk (all |a| < R), the identity (w − a) = b_R(a,w)·((R² − ā w)/R) gives
f = B·G with G(w) := g₀(w)·∏_{a∈S}((R² − ā w)/R)^{m a}:
* `lin_eq_b_mul`   — w − a = b R a w · ((R² − ā w)/R)   (|a| < R, |w| ≤ R);
* `fac_eq_bprod_mul` — (∏(w−a)^m)·g₀ w = B w · G w on the closed disk;
* `Gfac_ne_zero`   — g₀ w ≠ 0 ⇒ G w ≠ 0 (G is zero-free wherever g₀ is — ALL zeros of f in the disk sit in B);
* `Gfac_norm_sphere` — ‖w‖ = R ⇒ ‖G w‖ = ‖f w‖  (|B| = 1 on the sphere);
* `logDeriv_split` — f = B·G near 0, B 0 ≠ 0, G 0 ≠ 0, both differentiable at 0 ⇒ logDeriv f 0 = logDeriv B 0 + logDeriv G 0.
So sharp BC (S3018) can be run on G over the WHOLE radius R (not R/2), with sup‖G‖ = sup‖f‖ on the sphere.
-/

open Complex Finset Metric

namespace Blaschke

noncomputable def Gfac (R : ℝ) (S : Finset ℂ) (m : ℂ → ℕ) (g₀ : ℂ → ℂ) (w : ℂ) : ℂ :=
  g₀ w * ∏ a ∈ S, (((R : ℂ) ^ 2 - (starRingEnd ℂ) a * w) / (R : ℂ)) ^ m a

theorem lin_eq_b_mul {R : ℝ} {a w : ℂ} (ha : ‖a‖ < R) (hw : ‖w‖ ≤ R) :
    w - a = b R a w * (((R : ℂ) ^ 2 - (starRingEnd ℂ) a * w) / (R : ℂ)) := by
  have hR : 0 < R := lt_of_le_of_lt (norm_nonneg _) ha
  have hRc : (R : ℂ) ≠ 0 := by exact_mod_cast hR.ne'
  have hden := blaschke_den_ne ha hw
  have hden' : (R : ℂ) ^ 2 - w * (starRingEnd ℂ) a ≠ 0 := by rwa [mul_comm w]
  unfold b
  field_simp

theorem fac_eq_bprod_mul {R : ℝ} {S : Finset ℂ} (m : ℂ → ℕ) (g₀ : ℂ → ℂ) (hS : ∀ a ∈ S, ‖a‖ < R) {w : ℂ}
    (hw : ‖w‖ ≤ R) :
    (∏ a ∈ S, (w - a) ^ m a) * g₀ w = bprod R S m w * Gfac R S m g₀ w := by
  unfold bprod Gfac
  have h : ∏ a ∈ S, (w - a) ^ m a
      = (∏ a ∈ S, b R a w ^ m a) * ∏ a ∈ S, (((R : ℂ) ^ 2 - (starRingEnd ℂ) a * w) / (R : ℂ)) ^ m a := by
    rw [← Finset.prod_mul_distrib]
    refine Finset.prod_congr rfl (fun a ha => ?_)
    rw [← mul_pow, ← lin_eq_b_mul (hS a ha) hw]
  rw [h]; ring

theorem Gfac_ne_zero {R : ℝ} {S : Finset ℂ} (m : ℂ → ℕ) {g₀ : ℂ → ℂ} (hS : ∀ a ∈ S, ‖a‖ < R) {w : ℂ}
    (hw : ‖w‖ ≤ R) (hg : g₀ w ≠ 0) : Gfac R S m g₀ w ≠ 0 := by
  unfold Gfac
  refine mul_ne_zero hg (Finset.prod_ne_zero_iff.mpr (fun a ha => pow_ne_zero _ ?_))
  have hR : 0 < R := lt_of_le_of_lt (norm_nonneg _) (hS a ha)
  exact div_ne_zero (blaschke_den_ne (hS a ha) hw) (by exact_mod_cast hR.ne')

theorem Gfac_norm_sphere {R : ℝ} {S : Finset ℂ} (m : ℂ → ℕ) (g₀ : ℂ → ℂ) (hS : ∀ a ∈ S, ‖a‖ < R) {w : ℂ}
    (hw : ‖w‖ = R) :
    ‖Gfac R S m g₀ w‖ = ‖(∏ a ∈ S, (w - a) ^ m a) * g₀ w‖ := by
  rw [fac_eq_bprod_mul m g₀ hS hw.le, norm_mul, bprod_norm_sphere m hS hw, one_mul]

theorem logDeriv_split {f B G : ℂ → ℂ} (hfe : f =ᶠ[nhds 0] fun w => B w * G w)
    (hB : B 0 ≠ 0) (hG : G 0 ≠ 0) (hBd : DifferentiableAt ℂ B 0) (hGd : DifferentiableAt ℂ G 0) :
    logDeriv f 0 = logDeriv B 0 + logDeriv G 0 := by
  have h1 : logDeriv f 0 = logDeriv (fun w => B w * G w) 0 := by
    rw [logDeriv_apply, logDeriv_apply, hfe.deriv_eq, hfe.eq_of_nhds]
  rw [h1]
  exact logDeriv_mul 0 hB hG hBd hGd

end Blaschke

#print axioms Blaschke.lin_eq_b_mul
#print axioms Blaschke.fac_eq_bprod_mul
#print axioms Blaschke.Gfac_ne_zero
#print axioms Blaschke.Gfac_norm_sphere
#print axioms Blaschke.logDeriv_split
