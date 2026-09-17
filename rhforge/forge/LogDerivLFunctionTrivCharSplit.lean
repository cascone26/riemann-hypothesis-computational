import Mathlib.NumberTheory.LSeries.DirichletContinuation
import Mathlib.NumberTheory.LSeries.Dirichlet
import Mathlib.Analysis.Calculus.LogDeriv
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv

set_option maxHeartbeats 1000000

/-!
# §1553 — logDeriv split of the trivial-character L-function: `L(1) = (∏)·ζ` near `Re > 1`

Foundation for the QUADRATIC (`χ²=1`) case of the Dirichlet-L zero-free region.  When `χ²=1`,
`LFunction (χ²) = LFunction 1 = LFunctionTrivChar N`, which (Mathlib `LFunctionTrivChar_eq_mul_riemannZeta`) equals
`(∏_{p|N}(1-p^{-s}))·ζ(s)` for `s ≠ 1`.  On `Re s > 1` (so `s ≠ 1`, both factors nonzero and differentiable), the
functions agree on a neighborhood, so

  `logDeriv (LFunctionTrivChar N) s = logDeriv (∏_{p|N}(1-p^{-·})) s + logDeriv ζ s`.

The `ζ` term carries the 2iγ growth (handled by the ζ tower's own pole-aware machinery §1404); the finite product
term is a bounded correction.  EXIT0 / std-axioms.
-/

open Complex Filter DirichletCharacter
open scoped Topology

theorem logDeriv_LFunctionTrivChar_split {N : ℕ} [NeZero N] {s : ℂ} (hs : 1 < s.re) :
    logDeriv (LFunctionTrivChar N) s
      = logDeriv (fun w => ∏ p ∈ N.primeFactors, (1 - (p : ℂ) ^ (-w))) s
        + logDeriv riemannZeta s := by
  have hne : s ≠ 1 := by
    intro h; rw [h] at hs; simp at hs
  -- agreement on a neighborhood of s
  have heqOn : Set.EqOn (LFunctionTrivChar N)
      (fun w => (∏ p ∈ N.primeFactors, (1 - (p : ℂ) ^ (-w))) * riemannZeta w) {w : ℂ | w ≠ 1} :=
    fun w hw => LFunctionTrivChar_eq_mul_riemannZeta hw
  have hmem : {w : ℂ | w ≠ 1} ∈ 𝓝 s := isOpen_ne.mem_nhds hne
  have heq : LFunctionTrivChar N
      =ᶠ[𝓝 s] (fun w => (∏ p ∈ N.primeFactors, (1 - (p : ℂ) ^ (-w))) * riemannZeta w) :=
    heqOn.eventuallyEq_of_mem hmem
  rw [(logDeriv_congr_nhds heq).eq_of_nhds]
  -- product ≠ 0 at s
  have hprod_ne : (∏ p ∈ N.primeFactors, (1 - (p : ℂ) ^ (-s))) ≠ 0 := by
    rw [Finset.prod_ne_zero_iff]
    intro p hp
    have hp2 : 2 ≤ p := (Nat.prime_of_mem_primeFactors hp).two_le
    have hp0 : 0 < p := by omega
    have hnorm : ‖(p : ℂ) ^ (-s)‖ = (p : ℝ) ^ (-s.re) := by
      rw [Complex.norm_natCast_cpow_of_pos hp0, Complex.neg_re]
    have hlt : ‖(p : ℂ) ^ (-s)‖ < 1 := by
      rw [hnorm]
      apply Real.rpow_lt_one_of_one_lt_of_neg
      · exact_mod_cast hp2.trans_lt' (by norm_num)
      · simpa using lt_trans (by norm_num : (0:ℝ) < 1) hs
    intro hzero
    have : (p : ℂ) ^ (-s) = 1 := by linear_combination -hzero
    rw [this, norm_one] at hlt
    exact lt_irrefl 1 hlt
  have hζ_ne : riemannZeta s ≠ 0 := riemannZeta_ne_zero_of_one_lt_re hs
  have hdf : DifferentiableAt ℂ (fun w => ∏ p ∈ N.primeFactors, (1 - (p : ℂ) ^ (-w))) s := by
    have hrw : (fun w => ∏ p ∈ N.primeFactors, (1 - (p : ℂ) ^ (-w)))
        = ∏ p ∈ N.primeFactors, (fun w : ℂ => (1 - (p : ℂ) ^ (-w))) := by
      funext w; rw [Finset.prod_apply]
    rw [hrw]
    apply DifferentiableAt.finsetProd
    intro p hp
    have hp0 : (p : ℂ) ≠ 0 := by
      have := (Nat.prime_of_mem_primeFactors hp).two_le
      exact_mod_cast (by omega : p ≠ 0)
    have hexp : DifferentiableAt ℂ (fun w : ℂ => -w) s := differentiableAt_id.neg
    exact (differentiableAt_const 1).sub (hexp.const_cpow (Or.inl hp0))
  have hdg : DifferentiableAt ℂ riemannZeta s := differentiableAt_riemannZeta hne
  exact logDeriv_mul s hprod_ne hζ_ne hdf hdg

#print axioms logDeriv_LFunctionTrivChar_split
