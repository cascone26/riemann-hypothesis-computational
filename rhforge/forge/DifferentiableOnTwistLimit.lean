import UniformCauchySeqOnTwist
import TwistLimit
import DifferentiableTwistPartialSum
import Mathlib.Analysis.Calculus.UniformLimitsDeriv
import Mathlib.Analysis.Complex.CauchyIntegral

open Finset Filter Topology

/-! §1531 twistLimit χ is DifferentiableOn {Re>0} (holomorphic limit). EXIT0 / std-axioms. -/

theorem differentiableOn_twistLimit {M : ℕ} [NeZero M] {χ : DirichletCharacter ℂ M} (hχ : χ ≠ 1) :
    DifferentiableOn ℂ (twistLimit χ) {s : ℂ | 0 < s.re} := by
  have hU : IsOpen {s : ℂ | 0 < s.re} := isOpen_lt continuous_const Complex.continuous_re
  have hlu : TendstoLocallyUniformlyOn
      (fun n s => ∑ k ∈ Finset.Icc (1 : ℕ) n, χ (k : ZMod M) * (k : ℂ) ^ (-s))
      (twistLimit χ) atTop {s : ℂ | 0 < s.re} := by
    rw [tendstoLocallyUniformlyOn_iff_forall_isCompact hU]
    intro K hKsub hKcompact
    rcases K.eq_empty_or_nonempty with hKe | hKne
    · subst hKe; exact tendstoUniformlyOn_empty
    · obtain ⟨z₀, hz₀K, hz₀min⟩ :=
        hKcompact.exists_isMinOn hKne Complex.continuous_re.continuousOn
      obtain ⟨R, hR⟩ := hKcompact.exists_bound_of_continuousOn continuous_norm.continuousOn
      have hδ : 0 < z₀.re := hKsub hz₀K
      have hKreg : K ⊆ {s : ℂ | z₀.re ≤ s.re ∧ ‖s‖ ≤ R} := by
        intro s hs
        refine ⟨hz₀min hs, ?_⟩
        have := hR s hs
        rwa [Real.norm_of_nonneg (norm_nonneg s)] at this
      have hUC := (uniformCauchySeqOn_twist_partial_sum hχ hδ (R := R)).mono hKreg
      exact hUC.tendstoUniformlyOn_of_tendsto (fun s hs => tendsto_twistLimit hχ (hKsub hs))
  exact hlu.differentiableOn
    (Filter.Eventually.of_forall
      (fun n => (differentiable_twist_partial_sum χ n).differentiableOn)) hU

#print axioms differentiableOn_twistLimit
