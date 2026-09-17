import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.Analysis.Analytic.Constructions
import Mathlib.Analysis.Meromorphic.Divisor

/-!
# §1389 — the polynomial product `∏_{u∈T} (s-u)^{m u}` is analytic

The bare polynomial factor `P_half`/`P_rest` (no `g`): a finite product of `(s-u)^{m u}` is `AnalyticOnNhd` on any
ball.  Combined with §1383 (`P_rest · g` analytic), this gives the analyticity — hence continuity — of the full
split function `P_half · (P_rest · g)`, the last ingredient for the `hfac` sphere-pointwise upgrade §1385.
EXIT0 / std-axioms.
-/

open Metric

theorem finset_prod_sub_analyticOnNhd {T : Finset ℂ} {c : ℂ} {R : ℝ} (m : ℂ → ℕ) :
    AnalyticOnNhd ℂ (fun s => ∏ u ∈ T, (s - u) ^ (m u)) (ball c R) := by
  have heq : (fun s : ℂ => ∏ u ∈ T, (s - u) ^ (m u))
      = ∏ u ∈ T, (fun s : ℂ => (s - u) ^ (m u)) := by
    ext s; simp [Finset.prod_apply]
  rw [heq]
  apply Finset.analyticOnNhd_prod
  intro u _
  exact fun z _ => (analyticAt_id.sub analyticAt_const).pow (m u)

#print axioms finset_prod_sub_analyticOnNhd
