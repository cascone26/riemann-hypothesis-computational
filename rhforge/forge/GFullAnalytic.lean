import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.Analysis.Analytic.Constructions
import Mathlib.Analysis.Meromorphic.Divisor

/-!
# §1383 — `g_full = P_rest · g` is analytic on the disk

The polynomial `P_rest(s) = ∏_{ρ∈T} (s-ρ)^{m ρ}` is entire, and `g` is `AnalyticOnNhd` on `ball c R` (§1343), so
their product `g_full` is `AnalyticOnNhd` on `ball c R`.  Together with §1382 (nonvanishing on the inner disk)
this gives the analytic nonvanishing factor §1378 consumes.  EXIT0 / std-axioms.
-/

open Metric

theorem prod_rest_mul_analyticOnNhd {T : Finset ℂ} {c : ℂ} {R : ℝ} (m : ℂ → ℕ) {g : ℂ → ℂ}
    (hg : AnalyticOnNhd ℂ g (ball c R)) :
    AnalyticOnNhd ℂ (fun s => (∏ ρ ∈ T, (s - ρ) ^ (m ρ)) * g s) (ball c R) := by
  apply AnalyticOnNhd.mul _ hg
  have heq : (fun x : ℂ => ∏ ρ ∈ T, (x - ρ) ^ (m ρ))
      = ∏ ρ ∈ T, (fun x : ℂ => (x - ρ) ^ (m ρ)) := by
    ext x; simp [Finset.prod_apply]
  rw [heq]
  apply Finset.analyticOnNhd_prod
  intro ρ _
  exact fun z _ => (analyticAt_id.sub analyticAt_const).pow (m ρ)

#print axioms prod_rest_mul_analyticOnNhd
