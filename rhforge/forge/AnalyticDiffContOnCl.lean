import Mathlib.Analysis.Analytic.Constructions
import Mathlib.Analysis.Calculus.DiffContOnCl
import Mathlib.Analysis.SpecialFunctions.Complex.Circle

/-!
# §1384 — `AnalyticOnNhd` on a larger ball ⟹ `DiffContOnCl` on the smaller ball

The standard radius trick supplying §1378's `DiffContOnCl g (ball c R)` hypothesis.  If `f` is analytic on the
larger ball `ball c R'` (`R < R'`), then on `closedBall c R ⊆ ball c R'` it is analytic — hence differentiable on
`ball c R` and continuous on `closure (ball c R) = closedBall c R`, i.e. `DiffContOnCl ℂ f (ball c R)`.  Apply
§1343's factorization at radius `R'` to obtain the factor on a ball just larger than the working disk.
EXIT0 / std-axioms.
-/

open Metric

theorem analyticOnNhd_diffContOnCl {c : ℂ} {R R' : ℝ} (hR : 0 < R) (hRR' : R < R') {f : ℂ → ℂ}
    (hf : AnalyticOnNhd ℂ f (ball c R')) : DiffContOnCl ℂ f (ball c R) := by
  refine ⟨(hf.mono (ball_subset_ball hRR'.le)).differentiableOn, ?_⟩
  rw [closure_ball c hR.ne']
  exact (hf.mono (closedBall_subset_ball hRR')).continuousOn

#print axioms analyticOnNhd_diffContOnCl
