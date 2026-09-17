import RegionLFunctionR14Clean
import LogSharpCloseLFunction
import PoleConstBound

set_option maxHeartbeats 1000000

/-!
# §1551 ★★★★★ — SHARP unconditional Dirichlet-L zero-free region `β ≤ 1 - c/log|γ|`

L-analogue of ζ §1462.  Same tower as §1549 but with the OPTIMAL constant choice `B = C = 64·log|γ|`
(so `K = 3A₀+320·log|γ| ≍ log|γ|`), closing the two log bounds via §1550 (log-log absorption).  Result: for every
nontrivial character `χ` with `χ²≠1`, an absolute `A₀ ≥ 0` with, for every zero `β+iγ` of `LFunction χ` past an
explicit (modulus-dependent) threshold,

  `β ≤ 1 - 1/(14·(3A₀ + 4·(64 log|γ|) + 64 log|γ|))`   (`= 1 - c/log|γ|`, the classical de la Vallée Poussin region).

EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem LFunction_zero_free_region_sharp {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) (hχ2 : χ ^ 2 ≠ 1) :
    ∃ A₀ : ℝ, 0 ≤ A₀ ∧ ∀ γ β : ℝ, Real.exp 16 ≤ |γ| →
      2 * Real.log ((N : ℝ) * (4 + 6) * (3 + 6 * A₀ + 640)) ≤ Real.log |γ| →
      LFunction χ ((β : ℂ) + I * (γ : ℝ)) = 0 →
      β ≤ 1 - 1 / (14 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|)) := by
  obtain ⟨A₀, hA₀nn, hA₀bd⟩ := exists_logDeriv_LFunctionTrivChar₁_bound
  have hN1 : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast Nat.one_le_iff_ne_zero.mpr (NeZero.ne N)
  refine ⟨A₀, hA₀nn, ?_⟩
  intro γ β hγe hγthr hL
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
  -- e=2 threshold from given e=4 threshold (8 ≤ 10, log monotone)
  have hXpos : (0 : ℝ) < 3 + 6 * A₀ + 640 := by linarith
  have hbase_mono : (N : ℝ) * (2 + 6) * (3 + 6 * A₀ + 640) ≤ (N : ℝ) * (4 + 6) * (3 + 6 * A₀ + 640) := by
    nlinarith [hN1, hXpos, mul_nonneg (show (0:ℝ) ≤ (N:ℝ) by linarith) hXpos.le]
  have hlogmono : Real.log ((N : ℝ) * (2 + 6) * (3 + 6 * A₀ + 640))
      ≤ Real.log ((N : ℝ) * (4 + 6) * (3 + 6 * A₀ + 640)) :=
    Real.log_le_log (by positivity) hbase_mono
  have hγthr2 : 2 * Real.log ((N : ℝ) * (2 + 6) * (3 + 6 * A₀ + 640)) ≤ Real.log |γ| := by
    linarith [hlogmono, hγthr]
  -- the two sharp log bounds via §1550
  have hBclean : 32 * Real.log ((N : ℝ) * (2 * |γ| + 6)
      * (3 + 2 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|))) ≤ 64 * Real.log |γ| := by
    have h := log_sharp_close_LFunction (N' := (N : ℝ)) (a := A₀) (x := |γ|) (e := 2)
      hN1 hA₀nn (by norm_num) hγe hγthr2
    have heq : (3 : ℝ) + 2 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|)
        = 3 + 6 * A₀ + 640 * Real.log |γ| := by ring
    rw [heq]; exact h
  have hCclean : 32 * Real.log ((N : ℝ) * (2 * |2 * γ| + 6)
      * (3 + 2 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|))) ≤ 64 * Real.log |γ| := by
    have h := log_sharp_close_LFunction (N' := (N : ℝ)) (a := A₀) (x := |γ|) (e := 4)
      hN1 hA₀nn (by norm_num) hγe hγthr
    have heq : (N : ℝ) * (2 * |2 * γ| + 6)
          * (3 + 2 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|))
        = (N : ℝ) * (4 * |γ| + 6) * (3 + 6 * A₀ + 640 * Real.log |γ|) := by
      rw [abs_mul]; simp only [abs_two]; ring
    rw [heq]; exact h
  exact region_LFunction_R14_clean hχ hχ2 hγ hKbig hA₀nn hB hB hL hBclean hA_pole hCclean

#print axioms LFunction_zero_free_region_sharp
