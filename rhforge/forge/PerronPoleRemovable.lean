import PerronIntegrandResidue
import PerronIntegrandDiff
import Mathlib.NumberTheory.LSeries.Dirichlet
import Mathlib.Analysis.Complex.RemovableSingularity

set_option maxHeartbeats 1000000

/-!
# §1613 — the Perron pole is removable: `(s-1)·F(s)` extends analytically across `s=1` with value `x²/2`

`F(s) = -ζ'/ζ(s)·x^{s+1}/(s(s+1))` has a simple pole at `s=1`; multiplying by `(s-1)` removes it.  Concretely
`Ĥ := update (fun s => (s-1)·F(s)) 1 (x²/2)` is analytic at `1`: it is differentiable on a punctured neighbourhood
(§1607, off the singular loci — all avoided near `1` except the pole) and continuous at `1` because §1608 gives
`(s-1)·F(s) → x²/2 = Ĥ(1)`.  Riemann's removable-singularity theorem
(`analyticAt_of_differentiable_on_punctured_nhds_of_continuousAt`) then yields analyticity.  This packages the
simple pole with its residue `x²/2` as a genuine holomorphic germ — the object the contour shift transports across
`s=1` to deposit the PNT main term.  EXIT0 / std-axioms.
-/

open Complex Filter Topology

theorem analyticAt_sub_one_mul_perron_integrand {x : ℝ} (hx : 0 < x) :
    AnalyticAt ℂ
      (Function.update
        (fun s => (s - 1) * (-logDeriv riemannZeta s * (x : ℂ) ^ (s + 1) / (s * (s + 1))))
        1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1 := by
  set F : ℂ → ℂ := fun s => -logDeriv riemannZeta s * (x : ℂ) ^ (s + 1) / (s * (s + 1)) with hF
  set H : ℂ → ℂ := fun s => (s - 1) * F s with hH
  set Ĥ : ℂ → ℂ := Function.update H 1 ((x : ℂ) ^ (2 : ℂ) / 2) with hĤ
  -- eventual side conditions near 1
  have hprod_ne : ∀ᶠ s in 𝓝[≠] (1 : ℂ), (s - 1) * riemannZeta s ≠ 0 :=
    riemannZeta_residue_one.eventually_ne one_ne_zero
  have hζ_ev : ∀ᶠ s in 𝓝[≠] (1 : ℂ), riemannZeta s ≠ 0 := by
    filter_upwards [hprod_ne] with s hs; intro h; apply hs; rw [h, mul_zero]
  have hev0 : ∀ᶠ s in 𝓝[≠] (1 : ℂ), s ≠ 0 :=
    (isOpen_ne.eventually_mem (show (1 : ℂ) ≠ 0 by norm_num)).filter_mono nhdsWithin_le_nhds
  have hevm1 : ∀ᶠ s in 𝓝[≠] (1 : ℂ), s ≠ -1 :=
    (isOpen_ne.eventually_mem (show (1 : ℂ) ≠ -1 by norm_num)).filter_mono nhdsWithin_le_nhds
  -- differentiable on a punctured neighbourhood
  have hd : ∀ᶠ s in 𝓝[≠] (1 : ℂ), DifferentiableAt ℂ Ĥ s := by
    filter_upwards [hζ_ev, hev0, hevm1, self_mem_nhdsWithin] with s hζ hs0 hsm1 hmem
    have hs1 : s ≠ 1 := by simpa using hmem
    have heq : Ĥ =ᶠ[𝓝 s] H := by
      filter_upwards [isOpen_ne.mem_nhds hs1] with z hz
      exact Function.update_of_ne hz _ _
    rw [heq.differentiableAt_iff]
    have hFd : DifferentiableAt ℂ F s := differentiableAt_perron_integrand hx hs1 hζ hs0 hsm1
    exact (differentiableAt_id.sub_const 1).mul hFd
  -- continuous at 1 via the §1608 residue limit
  have hc : ContinuousAt Ĥ 1 := by
    rw [hĤ, continuousAt_update_same]
    exact perron_integrand_residue_one hx
  exact analyticAt_of_differentiable_on_punctured_nhds_of_continuousAt hd hc

#print axioms analyticAt_sub_one_mul_perron_integrand
