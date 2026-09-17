import LogDerivTrivChar1OneEq
import Mathlib.NumberTheory.LSeries.RiemannZeta
import Mathlib.Analysis.Complex.CauchyIntegral
import Mathlib.Analysis.Calculus.FDeriv.Analytic

set_option maxHeartbeats 1000000

/-!
# §1612 — the Euler–Stieltjes constant: `-ζ'/ζ(s) - 1/(s-1) → -logDeriv g(1)` as `s → 1`

The finite (constant) part of the Laurent expansion of `-ζ'/ζ` at its simple pole `s=1`.  From §1570,
`logDeriv g(s) = 1/(s-1) + logDeriv ζ(s)` with `g = (s-1)ζ` entire, so `-ζ'/ζ(s) - 1/(s-1) = -logDeriv g(s)`, and
`logDeriv g` is continuous at `1` (`g` entire, `g(1) ≠ 0`).  Hence the limit exists and equals `-logDeriv g(1)`.
This constant is exactly what a principal-part subtraction needs to make `F(s) - (x²/2)/(s-1)` removable at `s=1`,
letting the contour cross the pole while depositing the §1608 residue `x²/2`.  EXIT0 / std-axioms.
-/

open Complex Filter Topology DirichletCharacter

theorem neg_logDeriv_zeta_sub_pole_tendsto :
    Tendsto (fun s : ℂ => -logDeriv riemannZeta s - 1 / (s - 1)) (𝓝[≠] (1 : ℂ))
      (𝓝 (-logDeriv (LFunctionTrivChar₁ 1) 1)) := by
  -- `g = (s-1)ζ` entire with `g(1) ≠ 0` ⟹ `logDeriv g` continuous at 1
  have hg_diff : Differentiable ℂ (LFunctionTrivChar₁ 1) := differentiable_LFunctionTrivChar₁ 1
  have hg1_ne : LFunctionTrivChar₁ 1 1 ≠ 0 := LFunctionTrivChar₁_apply_one_ne_zero 1
  have hana : AnalyticOnNhd ℂ (LFunctionTrivChar₁ 1) Set.univ :=
    analyticOnNhd_univ_iff_differentiable.mpr hg_diff
  have hcont_dg : ContinuousAt (deriv (LFunctionTrivChar₁ 1)) 1 :=
    (hana.deriv 1 (Set.mem_univ 1)).continuousAt
  have hcont_g : ContinuousAt (LFunctionTrivChar₁ 1) 1 := hg_diff.continuous.continuousAt
  have hcont : ContinuousAt (fun s => -logDeriv (LFunctionTrivChar₁ 1) s) 1 := by
    have hd := (hcont_dg.div hcont_g hg1_ne).neg
    have hfun : (fun s => -logDeriv (LFunctionTrivChar₁ 1) s)
        = (fun s => -(deriv (LFunctionTrivChar₁ 1) s / LFunctionTrivChar₁ 1 s)) := by
      funext s; rw [logDeriv_apply]
    rw [hfun]; exact hd
  -- eventually `ζ ≠ 0` near 1
  have hprod_ne : ∀ᶠ s in 𝓝[≠] (1 : ℂ), (s - 1) * riemannZeta s ≠ 0 :=
    riemannZeta_residue_one.eventually_ne one_ne_zero
  have hζ_ev : ∀ᶠ s in 𝓝[≠] (1 : ℂ), riemannZeta s ≠ 0 := by
    filter_upwards [hprod_ne] with s hs
    intro h; apply hs; rw [h, mul_zero]
  -- the algebraic identity, eventually
  have hev : (fun s : ℂ => -logDeriv riemannZeta s - 1 / (s - 1))
      =ᶠ[𝓝[≠] (1 : ℂ)] (fun s => -logDeriv (LFunctionTrivChar₁ 1) s) := by
    filter_upwards [hζ_ev, self_mem_nhdsWithin] with s hζ hmem
    have hs1 : s ≠ 1 := by simpa using hmem
    have h70 := logDeriv_LFunctionTrivChar₁_one_eq hs1 hζ
    linear_combination h70
  exact (tendsto_congr' hev).mpr hcont.continuousWithinAt

#print axioms neg_logDeriv_zeta_sub_pole_tendsto
