import DiffNegLogDerivZeta
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv

set_option maxHeartbeats 1000000

/-!
# §1607 — holomorphy of the full smoothed-Perron integrand `F(s) = -ζ'/ζ(s)·x^{s+1}/(s(s+1))`

The integrand of the smoothed Perron formula (kernel `x^{s+1}/(s(s+1))` for `ψ₁`) is complex-differentiable at every
`s` avoiding the four singular loci: `s ≠ 1` and `ζ s ≠ 0` (the `-ζ'/ζ` factor, §1606), `s ≠ 0` and `s ≠ -1`
(the rational kernel).  `x^{s+1}` is entire in `s` for `x > 0` (`DifferentiableAt.const_cpow`, base `x ≠ 0`).  This
is the pointwise-holomorphy hypothesis for Mathlib's rectangle-contour vanishing lemma
`integral_boundary_rect_eq_zero_of_differentiable_on_off_countable`, applied (after subtracting the `s=1`
principal part) on a rectangle inside the zero-free region §1575.  EXIT0 / std-axioms.
-/

open Complex

theorem differentiableAt_perron_integrand {x : ℝ} (hx : 0 < x) {s : ℂ}
    (hs1 : s ≠ 1) (hζ : riemannZeta s ≠ 0) (hs0 : s ≠ 0) (hsm1 : s ≠ -1) :
    DifferentiableAt ℂ
      (fun z => -logDeriv riemannZeta z * (x : ℂ) ^ (z + 1) / (z * (z + 1))) s := by
  have hxne : (x : ℂ) ≠ 0 := by exact_mod_cast ne_of_gt hx
  have hlog : DifferentiableAt ℂ (fun z => -logDeriv riemannZeta z) s :=
    differentiableAt_neg_logDeriv_riemannZeta hs1 hζ
  have hcpow : DifferentiableAt ℂ (fun z : ℂ => (x : ℂ) ^ (z + 1)) s :=
    (differentiableAt_id.add_const 1).const_cpow (Or.inl hxne)
  have hnum : DifferentiableAt ℂ (fun z => -logDeriv riemannZeta z * (x : ℂ) ^ (z + 1)) s :=
    hlog.mul hcpow
  have hden : DifferentiableAt ℂ (fun z : ℂ => z * (z + 1)) s :=
    differentiableAt_id.mul (differentiableAt_id.add_const 1)
  have hs1ne : s + 1 ≠ 0 := fun h => hsm1 (eq_neg_of_add_eq_zero_left h)
  have hden_ne : s * (s + 1) ≠ 0 := mul_ne_zero hs0 hs1ne
  exact hnum.div hden hden_ne

#print axioms differentiableAt_perron_integrand
