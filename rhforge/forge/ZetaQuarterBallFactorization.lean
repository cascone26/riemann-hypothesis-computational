import ZetaFactorizationAtCenter
import ZetaDivisorFinite
import SmulFactorizedEval
import EqOfCodiscreteContinuous
import Mathlib.Analysis.Meromorphic.Divisor

set_option maxHeartbeats 2000000

open Complex Metric Filter

/-!
# §1681 ★★★ — THE ζ EVERYWHERE-FACTORIZATION in §1680's shape: `ζ = P·g` on the closed ball

The crux instantiation: on any ball `ball c R'` avoiding the pole `s = 1` (with `1 < Re c`),
`ζ` factors EVERYWHERE on `closedBall c R` (`R < R'`) as `(∏_{ρ∈S}(s-ρ)^{m ρ})·g(s)` where `S` =
the divisor support in the quarter-ball `closedBall c (R/4)` (all genuine ζ-ZEROS — the divisor is
nonnegative on pole-free balls), `m` = the multiplicities, and `g` (= annulus-zeros polynomial ×
the tower's nonvanishing factor `g₀`) is `DiffContOnCl` on `ball c R` and nonvanishing on
`ball c (R/4)`. Chain: the tower's §1343 codiscrete factorization + §1341
`eq_of_codiscreteWithin_of_continuousAt` (both sides continuous at every ball point once the
divisor's nonnegativity converts `zpow` to `pow`) + product filter-split. Everything §1680 needs.
-/

theorem exists_zeta_quarterBall_factorization {c : ℂ} {R R' : ℝ}
    (hR : 0 < R) (hRR' : R < R') (hc : 1 < c.re) (hball : (1 : ℂ) ∉ ball c R') :
    ∃ (S : Finset ℂ) (m : ℂ → ℕ) (g : ℂ → ℂ),
      (∀ ρ ∈ S, ρ ∈ closedBall c (R / 4)) ∧
      (∀ ρ ∈ S, riemannZeta ρ = 0) ∧
      (∀ s ∈ closedBall c R, riemannZeta s = (∏ ρ ∈ S, (s - ρ) ^ (m ρ)) * g s) ∧
      DiffContOnCl ℂ g (ball c R) ∧
      (∀ u ∈ ball c (R / 4), g u ≠ 0) := by
  classical
  have hR' : 0 < R' := hR.trans hRR'
  -- ζ is analytic on the pole-free ball; its divisor there is nonnegative
  have hζan : AnalyticOnNhd ℂ riemannZeta (ball c R') := by
    apply AnalyticOnNhd.mono analyticOn_riemannZeta
    intro z hz
    simp only [Set.mem_compl_iff, Set.mem_singleton_iff]
    intro h1
    exact hball (h1 ▸ hz)
  have hd0 : ∀ u : ℂ, 0 ≤ MeromorphicOn.divisor riemannZeta (ball c R') u := by
    intro u
    have := MeromorphicOn.AnalyticOnNhd.divisor_nonneg hζan u
    simpa using this
  set d : ℂ → ℤ := fun u => MeromorphicOn.divisor riemannZeta (ball c R') u with hddef
  have hdf := divisor_riemannZeta_ball_support_finite c R'
  obtain ⟨g₀, hga, hgn, hfac_cod⟩ := exists_factorization_riemannZeta_at_center hR' hc
  -- data
  refine ⟨hdf.toFinset.filter (· ∈ closedBall c (R / 4)),
    fun ρ => (d ρ).toNat,
    fun s => (∏ ρ ∈ hdf.toFinset.filter (fun u => ¬ u ∈ closedBall c (R / 4)),
      (s - ρ) ^ ((d ρ).toNat)) * g₀ s,
    ?_, ?_, ?_, ?_, ?_⟩
  · -- location
    intro ρ hρ
    exact (Finset.mem_filter.mp hρ).2
  · -- genuine zeros
    intro ρ hρ
    obtain ⟨hρsupp, hρquarter⟩ := Finset.mem_filter.mp hρ
    have hρball : ρ ∈ ball c R' := by
      rw [mem_closedBall] at hρquarter
      rw [mem_ball]
      linarith
    have hdne : d ρ ≠ 0 := by
      have := Set.Finite.mem_toFinset hdf |>.mp hρsupp
      rwa [Function.mem_support] at this
    have hdpos : 0 < d ρ := lt_of_le_of_ne (hd0 ρ) (Ne.symm hdne)
    have happ : d ρ = (meromorphicOrderAt riemannZeta ρ).untop₀ := by
      simp only [hddef]
      simp [MeromorphicOn.divisor_def, hζan.meromorphicOn, hρball]
    have horder : 0 < meromorphicOrderAt riemannZeta ρ := by
      by_contra hng
      push_neg at hng
      have hge : 0 ≤ meromorphicOrderAt riemannZeta ρ :=
        (hζan ρ hρball).meromorphicOrderAt_nonneg
      have heq0 : meromorphicOrderAt riemannZeta ρ = 0 := le_antisymm hng hge
      rw [happ, heq0] at hdpos
      simp at hdpos
    have htend := tendsto_zero_of_meromorphicOrderAt_pos horder
    have hcont : Tendsto riemannZeta (nhdsWithin ρ {ρ}ᶜ) (nhds (riemannZeta ρ)) :=
      ((hζan ρ hρball).continuousAt.tendsto).mono_left nhdsWithin_le_nhds
    exact tendsto_nhds_unique hcont htend
  · -- the everywhere-factorization
    intro s hscb
    have hsball' : s ∈ ball c R' := by
      rw [mem_closedBall] at hscb
      rw [mem_ball]
      linarith
    have hzpow : ∀ ρ : ℂ, ∀ s' : ℂ, (s' - ρ) ^ (d ρ) = (s' - ρ) ^ ((d ρ).toNat) := by
      intro ρ s'
      simp only [hddef]
      rw [← Int.toNat_of_nonneg (hd0 ρ), zpow_natCast, Int.toNat_of_nonneg (hd0 ρ)]
    have hFeq : ((∏ᶠ u, (· - u) ^ (d u)) • g₀ : ℂ → ℂ)
        = fun s' => (∏ ρ ∈ hdf.toFinset, (s' - ρ) ^ ((d ρ).toNat)) * g₀ s' := by
      funext s'
      rw [smul_factorized_rational_eval hdf s']
      congr 1
      exact Finset.prod_congr rfl (fun ρ _ => hzpow ρ s')
    have hζF : riemannZeta s
        = (∏ ρ ∈ hdf.toFinset, (s - ρ) ^ ((d ρ).toNat)) * g₀ s := by
      have h := eq_of_codiscreteWithin_of_continuousAt isOpen_ball hsball' hfac_cod
        (hζan s hsball').continuousAt ?_
      · rw [h, hFeq]
      · rw [hFeq]
        apply ContinuousAt.mul
        · apply Continuous.continuousAt
          apply continuous_finsetProd
          intro ρ _
          exact (continuous_id.sub continuous_const).pow _
        · exact (hga s hsball').continuousAt
    rw [hζF, ← Finset.prod_filter_mul_prod_filter_not hdf.toFinset (· ∈ closedBall c (R / 4))]
    ring
  · -- g is DiffContOnCl on the ball
    have hgd : DifferentiableOn ℂ
        (fun s => (∏ ρ ∈ hdf.toFinset.filter (fun u => ¬ u ∈ closedBall c (R / 4)),
          (s - ρ) ^ ((d ρ).toNat)) * g₀ s) (closedBall c R) := by
      apply DifferentiableOn.mul
      · apply Differentiable.differentiableOn
        apply Differentiable.fun_finsetProd
        intro ρ _
        exact (differentiable_id.sub_const ρ).pow _
      · intro s hs
        have hsball' : s ∈ ball c R' := by
          rw [mem_closedBall] at hs
          rw [mem_ball]
          linarith
        exact ((hga s hsball').differentiableAt).differentiableWithinAt
    rw [← closure_ball c hR.ne'] at hgd
    exact hgd.diffContOnCl
  · -- g nonvanishing on the quarter-ball
    intro u hu
    apply mul_ne_zero
    · apply Finset.prod_ne_zero_iff.mpr
      intro ρ hρ
      apply pow_ne_zero
      apply sub_ne_zero.mpr
      intro huρ
      apply (Finset.mem_filter.mp hρ).2
      rw [← huρ]
      exact ball_subset_closedBall hu
    · have huball : u ∈ ball c R' := by
        rw [mem_ball] at hu ⊢
        linarith
      exact hgn ⟨u, huball⟩

#print axioms exists_zeta_quarterBall_factorization
