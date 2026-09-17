import Mathlib.NumberTheory.LSeries.DirichletContinuation
import Mathlib.NumberTheory.LSeries.Dirichlet
import Mathlib.Analysis.Calculus.LogDeriv

set_option maxHeartbeats 1000000

/-!
# §1570 — `logDeriv (LFunctionTrivChar₁ 1) s = 1/(s-1) + logDeriv ζ s`  (`s ≠ 1`, `ζ s ≠ 0`)

Second brick of the SMALL-|t| arc.  `LFunctionTrivChar₁ 1 = (s-1)·ζ(s)` away from `1` (`Function.update_of_ne` +
`LFunctionTrivChar_eq_mul_riemannZeta` at `N=1`, `primeFactors 1 = ∅`), so its logDeriv splits by `logDeriv_mul`.
Rearranged: `-logDeriv ζ s = 1/(s-1) - logDeriv (LFunctionTrivChar₁ 1) s` — the COMPLEX-point pole decomposition
(§1405 gave only real `σ`).  With §1569 and a rectangle bound on `‖logDeriv (LFunctionTrivChar₁ 1)‖` this yields
the complex `hpole`:  `Re(-logDeriv ζ s) ≤ 1/(Re s - 1) + A₀` for `s` near `1` with `Re s > 1`, ANY `|Im s|`.
EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Filter
open scoped Topology

theorem logDeriv_LFunctionTrivChar₁_one_eq {s : ℂ} (hs1 : s ≠ 1) (hζ : riemannZeta s ≠ 0) :
    logDeriv (LFunctionTrivChar₁ 1) s = 1 / (s - 1) + logDeriv riemannZeta s := by
  have heqOn : Set.EqOn (LFunctionTrivChar₁ 1) (fun w => (w - 1) * riemannZeta w) {w : ℂ | w ≠ 1} := by
    intro w hw
    have hw' : w ≠ 1 := hw
    rw [LFunctionTrivChar₁, Function.update_of_ne hw']
    have : LFunctionTrivChar 1 w = riemannZeta w := by
      rw [LFunctionTrivChar_eq_mul_riemannZeta hw', Nat.primeFactors_one, Finset.prod_empty, one_mul]
    rw [this]
  have hmem : {w : ℂ | w ≠ 1} ∈ 𝓝 s := isOpen_ne.mem_nhds hs1
  have heq : LFunctionTrivChar₁ 1 =ᶠ[𝓝 s] (fun w => (w - 1) * riemannZeta w) :=
    heqOn.eventuallyEq_of_mem hmem
  rw [(logDeriv_congr_nhds heq).eq_of_nhds]
  have hsub : s - 1 ≠ 0 := sub_ne_zero.mpr hs1
  have hd1 : DifferentiableAt ℂ (fun w : ℂ => w - 1) s := differentiableAt_id.sub_const 1
  have hd2 : DifferentiableAt ℂ riemannZeta s := differentiableAt_riemannZeta hs1
  rw [logDeriv_mul s hsub hζ hd1 hd2]
  congr 1
  rw [logDeriv_apply]
  have hderiv : deriv (fun w : ℂ => w - 1) s = 1 := by
    rw [deriv_sub_const]; exact deriv_id s
  rw [hderiv]

#print axioms logDeriv_LFunctionTrivChar₁_one_eq
