import LogDerivTrivChar1OneEq
import Mathlib.NumberTheory.LSeries.RiemannZeta
import Mathlib.Analysis.Complex.CauchyIntegral
import Mathlib.Analysis.Calculus.FDeriv.Analytic

set_option maxHeartbeats 1000000

/-!
# §1603 — the residue of `-ζ'/ζ` at `s = 1` is `1`: `(s-1)·(-ζ'/ζ(s)) → 1` as `s → 1`

The analytic engine of the PNT main term.  From §1570, `-ζ'/ζ(s) = 1/(s-1) - logDeriv g(s)` with
`g = LFunctionTrivChar₁ 1 = (s-1)ζ(s)` entire and `g(1) ≠ 0`.  Hence
`(s-1)·(-ζ'/ζ(s)) = 1 - (s-1)·logDeriv g(s)`, and since `g` is entire with `g(1) ≠ 0`, `logDeriv g` is
continuous at `1`, so `(s-1)·logDeriv g(s) → 0` and the whole thing → `1`.  This is exactly the simple pole with
residue `1` of `-ζ'/ζ` at `s=1`; contour-shifting the smoothed Perron integral §1601 past it yields the `x²/2`
main term of `Ψ`.  EXIT0 / std-axioms.
-/

open Complex Filter Topology DirichletCharacter

theorem neg_logDeriv_zeta_residue_one :
    Tendsto (fun s : ℂ => (s - 1) * (-logDeriv riemannZeta s)) (𝓝[≠] (1 : ℂ)) (𝓝 1) := by
  -- `g = (s-1)ζ(s)` is entire with `g(1) ≠ 0`
  have hg_diff : Differentiable ℂ (LFunctionTrivChar₁ 1) := differentiable_LFunctionTrivChar₁ 1
  have hg1_ne : LFunctionTrivChar₁ 1 1 ≠ 0 := LFunctionTrivChar₁_apply_one_ne_zero 1
  have hana : AnalyticOnNhd ℂ (LFunctionTrivChar₁ 1) Set.univ :=
    analyticOnNhd_univ_iff_differentiable.mpr hg_diff
  have hcont_dg : ContinuousAt (deriv (LFunctionTrivChar₁ 1)) 1 :=
    (hana.deriv 1 (Set.mem_univ 1)).continuousAt
  have hcont_g : ContinuousAt (LFunctionTrivChar₁ 1) 1 := hg_diff.continuous.continuousAt
  -- `(s-1)·logDeriv g(s) → 0`
  have hlog_tendsto :
      Tendsto (fun s : ℂ => (s - 1) * logDeriv (LFunctionTrivChar₁ 1) s) (𝓝[≠] (1 : ℂ)) (𝓝 0) := by
    have hdiv : Tendsto (fun s : ℂ => deriv (LFunctionTrivChar₁ 1) s / LFunctionTrivChar₁ 1 s)
        (𝓝[≠] (1 : ℂ)) (𝓝 (deriv (LFunctionTrivChar₁ 1) 1 / LFunctionTrivChar₁ 1 1)) :=
      (hcont_dg.div hcont_g hg1_ne).continuousWithinAt
    have hs1_tendsto : Tendsto (fun s : ℂ => s - 1) (𝓝[≠] (1 : ℂ)) (𝓝 0) := by
      have : Tendsto (fun s : ℂ => s - 1) (𝓝[≠] (1 : ℂ)) (𝓝 (1 - 1)) :=
        (continuousAt_id.sub continuousAt_const).continuousWithinAt
      simpa using this
    have hmul := hs1_tendsto.mul hdiv
    rw [zero_mul] at hmul
    have hfun : (fun s : ℂ => (s - 1) * logDeriv (LFunctionTrivChar₁ 1) s)
        = (fun s : ℂ => (s - 1) * (deriv (LFunctionTrivChar₁ 1) s / LFunctionTrivChar₁ 1 s)) := by
      funext s; rw [logDeriv_apply]
    rw [hfun]; exact hmul
  -- eventually `ζ s ≠ 0` near `1` (residue is `1 ≠ 0`)
  have hprod_ne : ∀ᶠ s in 𝓝[≠] (1 : ℂ), (s - 1) * riemannZeta s ≠ 0 :=
    riemannZeta_residue_one.eventually_ne one_ne_zero
  have hζ_ev : ∀ᶠ s in 𝓝[≠] (1 : ℂ), riemannZeta s ≠ 0 := by
    filter_upwards [hprod_ne] with s hs
    intro h; apply hs; rw [h, mul_zero]
  -- the algebraic identity, eventually
  have heq_ev : (fun s : ℂ => (s - 1) * (-logDeriv riemannZeta s))
      =ᶠ[𝓝[≠] (1 : ℂ)] (fun s : ℂ => 1 - (s - 1) * logDeriv (LFunctionTrivChar₁ 1) s) := by
    filter_upwards [hζ_ev, self_mem_nhdsWithin] with s hζ hmem
    have hs1 : s ≠ 1 := by simpa using hmem
    have hsub : s - 1 ≠ 0 := sub_ne_zero.mpr hs1
    have h70 := logDeriv_LFunctionTrivChar₁_one_eq hs1 hζ
    have key : (s - 1) * logDeriv (LFunctionTrivChar₁ 1) s
        = 1 + (s - 1) * logDeriv riemannZeta s := by
      rw [h70, mul_add]
      congr 1
      rw [mul_one_div, div_self hsub]
    rw [key]; ring
  have hfin : Tendsto (fun s : ℂ => (1 : ℂ) - (s - 1) * logDeriv (LFunctionTrivChar₁ 1) s)
      (𝓝[≠] (1 : ℂ)) (𝓝 (1 - 0)) := tendsto_const_nhds.sub hlog_tendsto
  rw [sub_zero] at hfin
  exact (tendsto_congr' heq_ev).mpr hfin

#print axioms neg_logDeriv_zeta_residue_one
