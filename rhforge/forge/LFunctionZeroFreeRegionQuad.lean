import RegionLFunctionR14CleanQuad
import Log32CloseLFunction
import Log32CloseSlack
import SumLogPrimeFactorsLe
import PoleConstBound

set_option maxHeartbeats 1000000

/-!
# §1564 ★★★★★ — UNCONDITIONAL quantitative zero-free region for `LFunction χ` with `χ²=1`

Quadratic analogue of §1549.  With `A₀` the ζ pole constant (§1360) and `B = C = |γ|` (`K = 3A₀+5|γ|`): the χ side
`hBclean` via §1548 (e=2), the χ² side `hCclean` via §1563 (the slack estimate, `d=6`, `k = ∑_{p|N} log p`), where
the slack threshold uses `∑_{p|N} log p ≤ log N` (§1562).  Result: for a nontrivial `χ` with `χ²=1`,

  `β ≤ 1 - 1/(14·(3A₀ + 5|γ|))`   (`= 1 - c/|γ|`).

EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem LFunction_zero_free_region_quad {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) (hχ2eq : χ ^ 2 = 1) :
    ∃ A₀ : ℝ, 0 ≤ A₀ ∧ ∀ γ β : ℝ, 65536 ≤ |γ| →
      128 * Real.log ((N : ℝ) * (2 + 6) * (16 + 6 * A₀)) ≤ |γ| →
      128 * Real.log (16 + 6 * A₀) + 2 * Real.log N ≤ |γ| →
      LFunction χ ((β : ℂ) + I * (γ : ℝ)) = 0 →
      β ≤ 1 - 1 / (14 * (3 * A₀ + 4 * |γ| + |γ|)) := by
  obtain ⟨A₀, hA₀nn, hA₀bd⟩ := exists_logDeriv_LFunctionTrivChar₁_bound
  have hN1 : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast Nat.one_le_iff_ne_zero.mpr (NeZero.ne N)
  have hNpos : 1 ≤ N := Nat.one_le_iff_ne_zero.mpr (NeZero.ne N)
  refine ⟨A₀, hA₀nn, ?_⟩
  intro γ β hγ1 hγ2a hγ2b hL
  have hγ : (1 : ℝ) < |γ| := by linarith
  have hKbig : 32 / 7 ≤ 3 * A₀ + 4 * |γ| + |γ| := by nlinarith [hA₀nn, hγ1]
  have hK0 : 0 < 3 * A₀ + 4 * |γ| + |γ| := by linarith
  have hσKmem : 1 + 1 / (2 * (3 * A₀ + 4 * |γ| + |γ|)) ∈ Set.Icc (1 : ℝ) 2 := by
    constructor
    · have : (0 : ℝ) < 1 / (2 * (3 * A₀ + 4 * |γ| + |γ|)) := by positivity
      linarith
    · have : 1 / (2 * (3 * A₀ + 4 * |γ| + |γ|)) ≤ 1 := by
        rw [div_le_one (by positivity)]; linarith
      linarith
  have hA_pole : ‖logDeriv (LFunctionTrivChar₁ 1)
      ((((1 + 1 / (2 * (3 * A₀ + 4 * |γ| + |γ|)) : ℝ)) : ℂ))‖ ≤ A₀ :=
    hA₀bd _ hσKmem
  -- χ side clean bound via §1548 (e=2)
  have hBclean : 32 * Real.log ((N : ℝ) * (2 * |γ| + 6) * (3 + 2 * (3 * A₀ + 4 * |γ| + |γ|))) ≤ |γ| := by
    have h := log32_close_LFunction (N' := (N : ℝ)) (a := A₀) (x := |γ|) (e := 2)
      hN1 hA₀nn (by norm_num) hγ1 hγ2a
    have heq : (3 : ℝ) + 2 * (3 * A₀ + 4 * |γ| + |γ|) = 3 + 6 * A₀ + 10 * |γ| := by ring
    rw [heq]; exact h
  -- χ² side clean bound via §1563 (slack, d=6, k = ∑ log p)
  have hsum := sum_log_primeFactors_le_log hNpos
  have hCclean : 32 * Real.log ((3 * |2 * γ| + 10) * (3 + 2 * (3 * A₀ + 4 * |γ| + |γ|)))
      ≤ |γ| - ∑ p ∈ N.primeFactors, Real.log p := by
    have hthr : 128 * Real.log (16 + 6 * A₀) + 2 * (∑ p ∈ N.primeFactors, Real.log p) ≤ |γ| := by
      linarith [hsum, hγ2b]
    have h := log32_close_slack (a := A₀) (x := |γ|) (d := 6)
      (k := ∑ p ∈ N.primeFactors, Real.log p) hA₀nn (by norm_num) (by norm_num) hγ1 hthr
    have heq : (3 * |2 * γ| + 10) * (3 + 2 * (3 * A₀ + 4 * |γ| + |γ|))
        = (6 * |γ| + 10) * (3 + 6 * A₀ + 10 * |γ|) := by
      rw [abs_mul]; simp only [abs_two]; ring
    rw [heq]; exact h
  exact region_LFunction_R14_clean_quad hχ hχ2eq hγ hKbig hA₀nn (abs_nonneg γ) (abs_nonneg γ)
    hL hBclean hA_pole hCclean

#print axioms LFunction_zero_free_region_quad
