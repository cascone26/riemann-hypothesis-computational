import RegionLFunctionR14Clean
import Log32CloseLFunction
import PoleConstBound

set_option maxHeartbeats 1000000

/-!
# §1549 ★★★★★ — UNCONDITIONAL quantitative Dirichlet-L zero-free region

L-analogue of ζ §1460.  Assembles the whole L tower.  With `A₀` the absolute ζ pole constant (§1360) and the
choice `B = C = |γ|` (`K = 3A₀+5|γ|`), every hypothesis of §1547 is verified: `hA_pole` from §1360, and the two
log bounds `hBclean`/`hCclean` from §1548 (`e=2`, iγ; `e=4`, 2iγ, since `2·|2γ| = 4|γ|`).  Result: for every
nontrivial character `χ` with `χ²≠1`, there is an absolute constant `A₀ ≥ 0` such that for every zero `β+iγ` of
`LFunction χ` with `|γ|` beyond an explicit (modulus-dependent) threshold,

  `β ≤ 1 - 1/(14·(3A₀ + 5|γ|))`   (`= 1 - c/|γ|`, a genuine quantitative zero-free region).

EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem LFunction_zero_free_region {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) (hχ2 : χ ^ 2 ≠ 1) :
    ∃ A₀ : ℝ, 0 ≤ A₀ ∧ ∀ γ β : ℝ, 65536 ≤ |γ| →
      128 * Real.log ((N : ℝ) * (4 + 6) * (16 + 6 * A₀)) ≤ |γ| →
      LFunction χ ((β : ℂ) + I * (γ : ℝ)) = 0 →
      β ≤ 1 - 1 / (14 * (3 * A₀ + 4 * |γ| + |γ|)) := by
  obtain ⟨A₀, hA₀nn, hA₀bd⟩ := exists_logDeriv_LFunctionTrivChar₁_bound
  have hN1 : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast Nat.one_le_iff_ne_zero.mpr (NeZero.ne N)
  refine ⟨A₀, hA₀nn, ?_⟩
  intro γ β hγ1 hγ2 hL
  have hγ : (1 : ℝ) < |γ| := by linarith
  have hKbig : 32 / 7 ≤ 3 * A₀ + 4 * |γ| + |γ| := by nlinarith [hA₀nn, hγ1]
  have hK0 : 0 < 3 * A₀ + 4 * |γ| + |γ| := by linarith
  -- pole bound at σK
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
  -- e=4 threshold given; e=2 threshold follows (8 ≤ 10, log monotone)
  have hbase_mono : (N : ℝ) * (2 + 6) * (16 + 6 * A₀) ≤ (N : ℝ) * (4 + 6) * (16 + 6 * A₀) := by
    have hb : (0 : ℝ) < 16 + 6 * A₀ := by linarith
    nlinarith [hN1, hb]
  have hlogmono : Real.log ((N : ℝ) * (2 + 6) * (16 + 6 * A₀))
      ≤ Real.log ((N : ℝ) * (4 + 6) * (16 + 6 * A₀)) :=
    Real.log_le_log (by positivity) hbase_mono
  have hγ2' : 128 * Real.log ((N : ℝ) * (2 + 6) * (16 + 6 * A₀)) ≤ |γ| := by linarith [hlogmono, hγ2]
  -- the two clean log bounds via §1548
  have hBclean : 32 * Real.log ((N : ℝ) * (2 * |γ| + 6) * (3 + 2 * (3 * A₀ + 4 * |γ| + |γ|))) ≤ |γ| := by
    have h := log32_close_LFunction (N' := (N : ℝ)) (a := A₀) (x := |γ|) (e := 2)
      hN1 hA₀nn (by norm_num) hγ1 hγ2'
    have heq1 : (2 : ℝ) * |γ| + 6 = 2 * |γ| + 6 := rfl
    have heq2 : (3 : ℝ) + 2 * (3 * A₀ + 4 * |γ| + |γ|) = 3 + 6 * A₀ + 10 * |γ| := by ring
    rw [heq2]
    calc 32 * Real.log ((N : ℝ) * (2 * |γ| + 6) * (3 + 6 * A₀ + 10 * |γ|))
        = 32 * Real.log ((N : ℝ) * (2 * |γ| + 6) * (3 + 6 * A₀ + 10 * |γ|)) := rfl
      _ ≤ |γ| := h
  have hCclean : 32 * Real.log ((N : ℝ) * (2 * |2 * γ| + 6) * (3 + 2 * (3 * A₀ + 4 * |γ| + |γ|))) ≤ |γ| := by
    have h := log32_close_LFunction (N' := (N : ℝ)) (a := A₀) (x := |γ|) (e := 4)
      hN1 hA₀nn (by norm_num) hγ1 hγ2
    have heq : (N : ℝ) * (2 * |2 * γ| + 6) * (3 + 2 * (3 * A₀ + 4 * |γ| + |γ|))
        = (N : ℝ) * (4 * |γ| + 6) * (3 + 6 * A₀ + 10 * |γ|) := by
      rw [abs_mul]; simp only [abs_two]; ring
    rw [heq]; exact h
  exact region_LFunction_R14_clean hχ hχ2 hγ hKbig hA₀nn (abs_nonneg γ) (abs_nonneg γ)
    hL hBclean hA_pole hCclean

#print axioms LFunction_zero_free_region
