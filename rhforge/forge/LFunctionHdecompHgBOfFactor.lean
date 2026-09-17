import LFunctionHalfDiskDecomp
import HgBLFunctionConcreteZpow

/-!
# §1493 — coordinated `hdecomp ∧ hgB` from ONE §1476 factorization (shared `g` for the region composition)

L-analogue of ζ §1410 (cleaner: entire, no `hball`).  From a single §1476 factorization of `LFunction χ` on
`ball c R'` (plus the sphere bound `Mζ`, `R<R'`, and the strict factor bound `hgc_lt`), produces the **same**
factor value `g = logDeriv g_full c` serving BOTH analytic slots of the generic region composition §1336:
* `hdecomp` — the half-disk logarithmic-derivative decomposition (§1486, `p = (·∈closedBall c (R/2))`);
* `hgB` — the Borel–Carathéodory lower bound `-bound ≤ g.re` (§1488).

The two bounds are about the *identical* `g_full` built from the *identical* `g₀` — the coordination content.
`c` kept abstract (only `1 < c.re`).  This is the packaged input for the final L region assembly.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric Filter Classical

theorem LFunction_hdecomp_and_hgB_of_factor {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) {c : ℂ} {R' R Mζ : ℝ} (hR : 0 < R) (hRR' : R < R') (hcre : 1 < c.re)
    {g₀ : ℂ → ℂ} (hga : AnalyticOnNhd ℂ g₀ (ball c R')) (hgn : ∀ u : ball c R', g₀ u ≠ 0)
    (hfac_cod : LFunction χ =ᶠ[codiscreteWithin (ball c R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u)) • g₀)
    (hLsphere : ∀ s ∈ sphere c R, ‖LFunction χ s‖ ≤ Mζ)
    (hgc_lt : ‖(∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
        (fun u => ¬ u ∈ closedBall c (R / 2)),
        (c - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) * g₀ c‖
      < Mζ / (R / 2) ^ (∑ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
        (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat)) :
    ∃ g : ℂ,
      logDeriv (LFunction χ) c
        = (∑ ρ ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
            (fun u => u ∈ closedBall c (R / 2)),
            (MeromorphicOn.divisor (LFunction χ) (ball c R') ρ : ℂ) / (c - ρ)) + g
      ∧ -(4 * (Real.log (Mζ / (R / 2) ^ (∑ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
          (· ∈ closedBall c (R / 2)), (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat))
          - Real.log ‖(∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
            (fun u => ¬ u ∈ closedBall c (R / 2)),
            (c - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) * g₀ c‖) / (R / 2)) ≤ g.re := by
  have hR' : 0 < R' := lt_trans hR hRR'
  refine ⟨logDeriv (fun s => (∏ ρ ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter
      (fun u => ¬ u ∈ closedBall c (R / 2)),
      (s - ρ) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') ρ)) * g₀ s) c, ?_, ?_⟩
  · exact logDeriv_LFunction_halfdisk_decomp_of_factor hχ (mem_ball_self hR') hcre hga hgn hfac_cod
      (fun u => u ∈ closedBall c (R / 2))
  · exact hgB_LFunction_concrete_zpow hχ hR hRR' hga hgn hfac_cod hLsphere hgc_lt

#print axioms LFunction_hdecomp_and_hgB_of_factor
