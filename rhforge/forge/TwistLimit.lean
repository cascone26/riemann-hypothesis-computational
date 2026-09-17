import CauchySeqTwistPartialSum

open Finset Filter Topology

/-! §1528 the conditional twisted-series limit g'(s) and its pointwise convergence on {Re>0}. EXIT0 / std. -/

noncomputable def twistLimit {M : ℕ} [NeZero M] (χ : DirichletCharacter ℂ M) (s : ℂ) : ℂ :=
  limUnder atTop (fun n => ∑ k ∈ Finset.Icc (1 : ℕ) n, χ (k : ZMod M) * (k : ℂ) ^ (-s))

theorem tendsto_twistLimit {M : ℕ} [NeZero M] {χ : DirichletCharacter ℂ M} (hχ : χ ≠ 1)
    {s : ℂ} (hσ : 0 < s.re) :
    Tendsto (fun n => ∑ k ∈ Finset.Icc (1 : ℕ) n, χ (k : ZMod M) * (k : ℂ) ^ (-s)) atTop
      (𝓝 (twistLimit χ s)) :=
  (cauchySeq_twist_partial_sum hχ hσ).tendsto_limUnder

#print axioms tendsto_twistLimit
