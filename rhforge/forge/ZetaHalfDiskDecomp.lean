import ZetaDecompAt
import LogDerivZpowMulSplit
import SubNeZeroDivisor

/-!
# §1393 — half-disk logDeriv decomposition of ζ (corrected §1349-analogue)

The decomposition the two-radii region actually needs: keep only the zeros satisfying a predicate `p` (the
half-disk zeros) in the partial-fractions sum; the rest (annulus zeros) are folded into a factor `g_full` whose
`Re logDeriv` is bounded by the two-radii §1378.  Obtained by regrouping §1349's *full* decomposition
(`logDeriv ζ = ∑_all divisor/(z-ρ) + logDeriv g_349`): split `∑_all = ∑_p + ∑_{¬p}`, and §1332 gives
`logDeriv g_full = ∑_{¬p} divisor/(z-ρ) + logDeriv g_349` where `g_full = (∏_{¬p}(·-ρ)^{divisor})·g_349`.  Hence

  `logDeriv ζ z = ∑_{p} divisor/(z-ρ) + logDeriv g_full z`.

This is `§1336`'s `hdecomp` with `S = half-disk zeros`, `g = logDeriv g_full z`.  EXIT0 / std-axioms.
-/

open Complex Metric Classical

theorem logDeriv_riemannZeta_halfdisk_decomp {c : ℂ} {R : ℝ} {z z₀ : ℂ}
    (hz : z ∈ ball c R) (hzre : 1 < z.re) (hz₀ : z₀ ∈ ball c R) (hz₀1 : z₀ ≠ 1)
    (hz₀ne : riemannZeta z₀ ≠ 0) (p : ℂ → Prop) :
    ∃ g : ℂ → ℂ, AnalyticOnNhd ℂ g (ball c R) ∧ (∀ u : ball c R, g u ≠ 0) ∧
      logDeriv riemannZeta z =
        (∑ ρ ∈ (divisor_riemannZeta_ball_support_finite c R).toFinset.filter p,
          (MeromorphicOn.divisor riemannZeta (ball c R) ρ : ℂ) / (z - ρ))
        + logDeriv (fun s => (∏ ρ ∈ (divisor_riemannZeta_ball_support_finite c R).toFinset.filter
            (fun u => ¬ p u), (s - ρ) ^ (MeromorphicOn.divisor riemannZeta (ball c R) ρ))
          * g s) z := by
  obtain ⟨g, hga, hgn, hdec⟩ := logDeriv_riemannZeta_decomp_at hz hzre hz₀ hz₀1 hz₀ne
  set D := MeromorphicOn.divisor riemannZeta (ball c R) with hD
  set T := (divisor_riemannZeta_ball_support_finite c R).toFinset with hT
  -- sub-nonzero on the whole support (for §1332)
  have hsub : ∀ ρ ∈ T.filter (fun u => ¬ p u), z - ρ ≠ 0 := by
    intro ρ hρ
    rw [Finset.mem_filter, hT, Set.Finite.mem_toFinset, Function.mem_support] at hρ
    have hρU : ρ ∈ ball c R := D.supportWithinDomain hρ.1
    exact sub_ne_zero_of_divisor_ne_zero hρU hρ.1 hzre
  have hgz : g z ≠ 0 := hgn ⟨z, hz⟩
  have hgd : DifferentiableAt ℂ g z := (hga z hz).differentiableAt
  -- §1332: logDeriv g_full z = ∑_{¬p} D/(z-ρ) + logDeriv g z
  have hsplit := logDeriv_prod_zpow_mul_split (S := T.filter (fun u => ¬ p u)) (m := D)
    (h := g) (z := z) hsub hgz hgd
  refine ⟨g, hga, hgn, ?_⟩
  rw [hdec, ← Finset.sum_filter_add_sum_filter_not T p, hsplit]
  ring

#print axioms logDeriv_riemannZeta_halfdisk_decomp
