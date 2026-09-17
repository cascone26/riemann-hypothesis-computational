import FinsetProdSubAnalytic
import GFullAnalytic

/-!
# §1390 — the split function `P_half · (P_rest · g)` is analytic

`AnalyticOnNhd.mul` of §1389 (`P_half = ∏_{T1}(s-u)^{m u}` analytic) and §1383 (`P_rest · g` analytic) gives the
full split function analytic on `ball c R` — hence continuous at each interior point, the continuity §1385 needs to
upgrade §1388's codiscrete factorization to pointwise `ζ z = P_half(z)·(P_rest(z)·g z)` on the sphere.
EXIT0 / std-axioms.
-/

open Metric

theorem split_fun_analyticOnNhd {T1 T2 : Finset ℂ} {c : ℂ} {R : ℝ} (m : ℂ → ℕ) {g : ℂ → ℂ}
    (hg : AnalyticOnNhd ℂ g (ball c R)) :
    AnalyticOnNhd ℂ
      (fun s => (∏ u ∈ T1, (s - u) ^ (m u)) * ((∏ u ∈ T2, (s - u) ^ (m u)) * g s)) (ball c R) :=
  (finset_prod_sub_analyticOnNhd m).mul (prod_rest_mul_analyticOnNhd m hg)

#print axioms split_fun_analyticOnNhd
