import ZetaOrderNeTop
import ZetaDivisorFinite
import Mathlib.Analysis.Meromorphic.FactorizedRational
import Mathlib.Analysis.Meromorphic.Divisor

/-!
# §1343 — ζ divisor factorization on a disk: `ζ = (∏ᶠ divisor)·g`, `g` analytic nonvanishing

Applying Mathlib's `MeromorphicOn.extract_zeros_poles` to ζ on `ball c R`.  ζ is meromorphic there
(`meromorphicOn_riemannZeta`), has finite order everywhere (§1330, ζ not locally zero — needs one nonvanishing
point in the connected ball), and finite divisor support (§1331), so it factors as
`ζ =ᶠ[codiscreteWithin (ball c R)] (∏ᶠ u,(·-u)^{divisor ζ u})•g` with `g` analytic and nonvanishing on the ball.
This supplies the `h` hypothesis of §1342 — the divisor factorization of ζ, obtained entirely from Mathlib.
EXIT0 / std-axioms.
-/

open Complex Metric Filter

theorem exists_factorization_riemannZeta {c : ℂ} {R : ℝ} {z₀ : ℂ}
    (hz₀ : z₀ ∈ ball c R) (hz₀1 : z₀ ≠ 1) (hz₀ne : riemannZeta z₀ ≠ 0) :
    ∃ g : ℂ → ℂ, AnalyticOnNhd ℂ g (ball c R) ∧ (∀ u : ball c R, g u ≠ 0) ∧
      riemannZeta =ᶠ[codiscreteWithin (ball c R)]
        (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R) u)) • g := by
  have hconn : IsConnected (ball c R) := ⟨⟨z₀, hz₀⟩, (convex_ball c R).isPreconnected⟩
  have hmero : MeromorphicOn riemannZeta (ball c R) :=
    meromorphicOn_riemannZeta.mono_set (Set.subset_univ _)
  have h₂ : ∀ u : ball c R, meromorphicOrderAt riemannZeta u ≠ ⊤ :=
    meromorphicOrderAt_riemannZeta_ne_top hconn hz₀ hz₀1 hz₀ne
  have h₃ : (MeromorphicOn.divisor riemannZeta (ball c R)).support.Finite :=
    divisor_riemannZeta_ball_support_finite c R
  exact hmero.extract_zeros_poles h₂ h₃

#print axioms exists_factorization_riemannZeta
