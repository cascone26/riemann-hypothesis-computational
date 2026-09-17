import ZetaDecompAt
import LogDerivZpowMulSplit

/-!
# §1400 — half-disk logDeriv decomposition of ζ with **external** factorization `g`

The coordination brick: like §1393 but taking the §1343 factorization data `(g, hga, hgn, hfac)` as *hypotheses*
instead of obtaining it internally.  This lets the final assembly obtain §1343 **once** and share the same `g`
(hence the same `g_full`) across both the half-disk decomposition (this brick, feeding §1336's `hdecomp`) and the
two-radii BC bound (§1391/§1378, feeding `hgB`) — resolving the g-coordination the region needs.  Proof = §1349's
fact-discharge + §1342 (external factorization) + Finset split + §1332 (regroup annulus into `g_full`).
EXIT0 / std-axioms.
-/

open Complex Metric Filter Classical

theorem logDeriv_riemannZeta_halfdisk_decomp_of_factor {c : ℂ} {R : ℝ} {z : ℂ}
    (hz : z ∈ ball c R) (hzre : 1 < z.re) {g : ℂ → ℂ}
    (hga : AnalyticOnNhd ℂ g (ball c R)) (hgn : ∀ u : ball c R, g u ≠ 0)
    (hfac : riemannZeta =ᶠ[codiscreteWithin (ball c R)]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor riemannZeta (ball c R) u)) • g)
    (p : ℂ → Prop) :
    logDeriv riemannZeta z =
      (∑ ρ ∈ (divisor_riemannZeta_ball_support_finite c R).toFinset.filter p,
        (MeromorphicOn.divisor riemannZeta (ball c R) ρ : ℂ) / (z - ρ))
      + logDeriv (fun s => (∏ ρ ∈ (divisor_riemannZeta_ball_support_finite c R).toFinset.filter
          (fun u => ¬ p u), (s - ρ) ^ (MeromorphicOn.divisor riemannZeta (ball c R) ρ)) * g s) z := by
  set T := (divisor_riemannZeta_ball_support_finite c R).toFinset with hT
  set D := MeromorphicOn.divisor riemannZeta (ball c R) with hD
  have hz1 : z ≠ 1 := fun h => by rw [h] at hzre; simp at hzre
  have hζc : ContinuousAt (logDeriv riemannZeta) z :=
    continuousAt_logDeriv (analyticAt_riemannZeta hz1) (riemannZeta_ne_zero_of_one_lt_re hzre)
  have hdz : D z = 0 := divisor_riemannZeta_eq_zero_of_one_lt_re hz hzre
  have hga_z : AnalyticAt ℂ g z := hga z hz
  have hgz : g z ≠ 0 := hgn ⟨z, hz⟩
  obtain ⟨hFa, hFne⟩ := analyticAt_and_ne_zero_factor hdz hga_z hgz
  have hFc : ContinuousAt (logDeriv ((∏ᶠ u, (· - u) ^ (D u)) • g)) z :=
    continuousAt_logDeriv hFa hFne
  have hsub : ∀ ρ ∈ T, z - ρ ≠ 0 := by
    intro ρ hρ
    rw [hT, Set.Finite.mem_toFinset, Function.mem_support] at hρ
    exact sub_ne_zero_of_divisor_ne_zero (D.supportWithinDomain hρ) hρ hzre
  have hfull := logDeriv_eq_sum_of_codiscrete_factor isOpen_ball hz
    (divisor_riemannZeta_ball_support_finite c R) hfac hζc hFc hsub hgz hga_z.differentiableAt
  have hsub' : ∀ ρ ∈ T.filter (fun u => ¬ p u), z - ρ ≠ 0 :=
    fun ρ hρ => hsub ρ (Finset.mem_of_mem_filter ρ hρ)
  have hsplit := logDeriv_prod_zpow_mul_split (S := T.filter (fun u => ¬ p u)) (m := D)
    (h := g) (z := z) hsub' hgz hga_z.differentiableAt
  rw [hfull, ← Finset.sum_filter_add_sum_filter_not T p, hsplit]
  ring

#print axioms logDeriv_riemannZeta_halfdisk_decomp_of_factor
