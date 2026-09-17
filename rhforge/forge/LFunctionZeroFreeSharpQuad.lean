import RegionLFunctionR14CleanQuad
import LogSharpCloseLFunction
import LogSharpCloseSlack
import SumLogPrimeFactorsLe
import PoleConstBound

set_option maxHeartbeats 1000000

/-!
# §1566 ★★★★★ — SHARP unconditional zero-free region for `LFunction χ` with `χ²=1`: `β ≤ 1 - c/log|γ|`

Quadratic analogue of §1551.  With `B = C = 64·log|γ|` (`K = 3A₀+320·log|γ|`): the χ side `hBclean` via §1550
(e=2), the χ² side `hCclean` via §1565 (sharp slack, d=6, `k = ∑_{p|N} log p`), whose threshold uses
`∑_{p|N} log p ≤ log N` (§1562).  Result: for nontrivial `χ` with `χ²=1`, past explicit (modulus-dependent)
thresholds,

  `β ≤ 1 - 1/(14·(3A₀ + 320·log|γ|))`   (`= 1 - c/log|γ|`, the classical de la Vallée Poussin order).

EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem LFunction_zero_free_region_sharp_quad {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) (hχ2eq : χ ^ 2 = 1) :
    ∃ A₀ : ℝ, 0 ≤ A₀ ∧ ∀ γ β : ℝ, Real.exp 16 ≤ |γ| →
      2 * Real.log ((N : ℝ) * (2 + 6) * (3 + 6 * A₀ + 640)) ≤ Real.log |γ| →
      2 * Real.log (16 * (3 + 6 * A₀ + 640)) + Real.log N / 16 ≤ Real.log |γ| →
      LFunction χ ((β : ℂ) + I * (γ : ℝ)) = 0 →
      β ≤ 1 - 1 / (14 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|)) := by
  obtain ⟨A₀, hA₀nn, hA₀bd⟩ := exists_logDeriv_LFunctionTrivChar₁_bound
  have hN1 : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast Nat.one_le_iff_ne_zero.mpr (NeZero.ne N)
  have hNpos : 1 ≤ N := Nat.one_le_iff_ne_zero.mpr (NeZero.ne N)
  refine ⟨A₀, hA₀nn, ?_⟩
  intro γ β hγe hγthr hγthr2 hL
  have he1 : (1 : ℝ) < Real.exp 16 := by nlinarith [Real.add_one_le_exp (16 : ℝ)]
  have hγ : (1 : ℝ) < |γ| := by linarith
  have hL16 : 16 ≤ Real.log |γ| := by
    calc (16 : ℝ) = Real.log (Real.exp 16) := (Real.log_exp 16).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos 16) hγe
  have hLpos : 0 < Real.log |γ| := by linarith
  have hB : 0 ≤ 64 * Real.log |γ| := by linarith
  have hK0 : 0 < 3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ| := by nlinarith [hA₀nn, hLpos]
  have hKbig : 32 / 7 ≤ 3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ| := by
    nlinarith [hA₀nn, hL16]
  have hσKmem : 1 + 1 / (2 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|))
      ∈ Set.Icc (1 : ℝ) 2 := by
    constructor
    · have : (0 : ℝ) < 1 / (2 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|)) := by positivity
      linarith
    · have : 1 / (2 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|)) ≤ 1 := by
        rw [div_le_one (by positivity)]; linarith
      linarith
  have hA_pole : ‖logDeriv (LFunctionTrivChar₁ 1)
      ((((1 + 1 / (2 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|)) : ℝ)) : ℂ))‖ ≤ A₀ :=
    hA₀bd _ hσKmem
  have hsum := sum_log_primeFactors_le_log hNpos
  -- χ side sharp clean bound via §1550 (e=2)
  have hBclean : 32 * Real.log ((N : ℝ) * (2 * |γ| + 6)
      * (3 + 2 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|))) ≤ 64 * Real.log |γ| := by
    have h := log_sharp_close_LFunction (N' := (N : ℝ)) (a := A₀) (x := |γ|) (e := 2)
      hN1 hA₀nn (by norm_num) hγe hγthr
    have heq : (3 : ℝ) + 2 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|)
        = 3 + 6 * A₀ + 640 * Real.log |γ| := by ring
    rw [heq]; exact h
  -- χ² side sharp clean bound via §1565 (slack, d=6, k = ∑ log p)
  have hCclean : 32 * Real.log ((3 * |2 * γ| + 10)
      * (3 + 2 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|)))
      ≤ 64 * Real.log |γ| - ∑ p ∈ N.primeFactors, Real.log p := by
    have hdiv16 : (∑ p ∈ N.primeFactors, Real.log p) / 16 ≤ Real.log N / 16 := by linarith [hsum]
    have hthr' : 2 * Real.log (16 * (3 + 6 * A₀ + 640))
        + (∑ p ∈ N.primeFactors, Real.log p) / 16 ≤ Real.log |γ| := by
      linarith [hdiv16, hγthr2]
    have h := log_sharp_close_slack (a := A₀) (x := |γ|) (d := 6)
      (k := ∑ p ∈ N.primeFactors, Real.log p) hA₀nn (by norm_num) (by norm_num) hγe hthr'
    have heq : (3 * |2 * γ| + 10)
          * (3 + 2 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|))
        = (6 * |γ| + 10) * (3 + 6 * A₀ + 640 * Real.log |γ|) := by
      rw [abs_mul]; simp only [abs_two]; ring
    rw [heq]; exact h
  exact region_LFunction_R14_clean_quad hχ hχ2eq hγ hKbig hA₀nn hB hB hL hBclean hA_pole hCclean

#print axioms LFunction_zero_free_region_sharp_quad
