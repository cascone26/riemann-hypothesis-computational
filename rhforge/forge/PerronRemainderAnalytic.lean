import PerronPoleRemovable
import Mathlib.Analysis.Calculus.DSlope

set_option maxHeartbeats 1000000

/-!
# §1615 — the holomorphic remainder `R = dslope Ĥ 1` is analytic at `s = 1`

The remainder of the principal-part split (§1614) is `R = dslope Ĥ 1` with `Ĥ = update ((s-1)F) 1 (x²/2)`.  Since
`Ĥ` is analytic at `1` (§1613), `R` is continuous at `1` (`continuousAt_dslope_same`, needing `Ĥ` differentiable
there) and differentiable on a punctured neighbourhood (`differentiableAt_dslope_of_ne` + §1607), so Riemann's
removable-singularity theorem makes `R` analytic at `1`.  Unlike the raw integrand `F`, the remainder has NO pole at
`1` — so its rectangle contour integral vanishes even across `s=1`, and the entire residue `x²/2` comes solely from
the explicit pole term `(x²/2)/(s-1)`.  EXIT0 / std-axioms.
-/

open Complex Filter Topology

theorem analyticAt_perron_remainder {x : ℝ} (hx : 0 < x) :
    AnalyticAt ℂ
      (dslope
        (Function.update
          (fun z => (z - 1) * (-logDeriv riemannZeta z * (x : ℂ) ^ (z + 1) / (z * (z + 1))))
          1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1) 1 := by
  set Ĥ : ℂ → ℂ :=
    Function.update
      (fun z => (z - 1) * (-logDeriv riemannZeta z * (x : ℂ) ^ (z + 1) / (z * (z + 1))))
      1 ((x : ℂ) ^ (2 : ℂ) / 2) with hĤ
  have hĤdiff : DifferentiableAt ℂ Ĥ 1 := (analyticAt_sub_one_mul_perron_integrand hx).differentiableAt
  have hc : ContinuousAt (dslope Ĥ 1) 1 := continuousAt_dslope_same.mpr hĤdiff
  -- eventual side conditions near 1
  have hprod_ne : ∀ᶠ s in 𝓝[≠] (1 : ℂ), (s - 1) * riemannZeta s ≠ 0 :=
    riemannZeta_residue_one.eventually_ne one_ne_zero
  have hζ_ev : ∀ᶠ s in 𝓝[≠] (1 : ℂ), riemannZeta s ≠ 0 := by
    filter_upwards [hprod_ne] with s hs; intro h; apply hs; rw [h, mul_zero]
  have hev0 : ∀ᶠ s in 𝓝[≠] (1 : ℂ), s ≠ 0 :=
    (isOpen_ne.eventually_mem (show (1 : ℂ) ≠ 0 by norm_num)).filter_mono nhdsWithin_le_nhds
  have hevm1 : ∀ᶠ s in 𝓝[≠] (1 : ℂ), s ≠ -1 :=
    (isOpen_ne.eventually_mem (show (1 : ℂ) ≠ -1 by norm_num)).filter_mono nhdsWithin_le_nhds
  have hd : ∀ᶠ s in 𝓝[≠] (1 : ℂ), DifferentiableAt ℂ (dslope Ĥ 1) s := by
    filter_upwards [hζ_ev, hev0, hevm1, self_mem_nhdsWithin] with s hζ hs0 hsm1 hmem
    have hs1 : s ≠ 1 := by simpa using hmem
    rw [differentiableAt_dslope_of_ne hs1]
    have heq : Ĥ =ᶠ[𝓝 s]
        (fun z => (z - 1) * (-logDeriv riemannZeta z * (x : ℂ) ^ (z + 1) / (z * (z + 1)))) := by
      filter_upwards [isOpen_ne.mem_nhds hs1] with z hz; rw [hĤ]; exact Function.update_of_ne hz _ _
    rw [heq.differentiableAt_iff]
    exact (differentiableAt_id.sub_const 1).mul (differentiableAt_perron_integrand hx hs1 hζ hs0 hsm1)
  exact analyticAt_of_differentiable_on_punctured_nhds_of_continuousAt hd hc

#print axioms analyticAt_perron_remainder
