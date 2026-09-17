import RegionR14Clean
import Log32Close
import BetaLtSigmaK
import PoleConstBound

/-!
# §1460 ★★★★★ — UNCONDITIONAL quantitative ζ zero-free region

Assembles the whole tower.  With `A₀` the absolute pole constant (§1360) and the choice `B = C = |γ|`
(`K = 3A₀+5|γ|`), every hypothesis of §1457 is verified: `hβ` from §1455, `hA_pole` from §1360, and the two
log bounds `hBclean`/`hCclean` from §1459 (`d=3`, `d=6`).  Result: there is an absolute constant `A₀ ≥ 0` such that
for every ζ zero `β+iγ` with `|γ|` beyond an explicit threshold,

  `β ≤ 1 - 1/(14·(3A₀ + 5|γ|))`   (`= 1 - c/|γ|`, a genuine quantitative zero-free region).

EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem riemannZeta_zero_free_region :
    ∃ A₀ : ℝ, 0 ≤ A₀ ∧ ∀ γ β : ℝ, 65536 ≤ |γ| → 128 * Real.log (16 + 6 * A₀) ≤ |γ| →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 →
      β ≤ 1 - 1 / (14 * (3 * A₀ + 4 * |γ| + |γ|)) := by
  obtain ⟨A₀, hA₀nn, hA₀bd⟩ := exists_logDeriv_LFunctionTrivChar₁_bound
  refine ⟨A₀, hA₀nn, ?_⟩
  intro γ β hγ1 hγ2 hζ
  have hγ : (1 : ℝ) < |γ| := by linarith
  have hKbig : 32 / 7 ≤ 3 * A₀ + 4 * |γ| + |γ| := by nlinarith [hA₀nn, hγ1]
  have hK0 : 0 < 3 * A₀ + 4 * |γ| + |γ| := by linarith
  -- pole bound at σK = 1+1/(2K)
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
  -- the two log bounds via §1459
  have hBclean : 32 * Real.log ((3 * |γ| + 10) * (3 + 2 * (3 * A₀ + 4 * |γ| + |γ|))) ≤ |γ| := by
    have h := log32_close (a := A₀) (x := |γ|) (d := 3) hA₀nn (by norm_num) (by norm_num) hγ1 hγ2
    have heq : (3 + 2 * (3 * A₀ + 4 * |γ| + |γ|)) = 3 + 6 * A₀ + 10 * |γ| := by ring
    rw [heq]; exact h
  have hCclean : 32 * Real.log ((3 * |2 * γ| + 10) * (3 + 2 * (3 * A₀ + 4 * |γ| + |γ|))) ≤ |γ| := by
    have h := log32_close (a := A₀) (x := |γ|) (d := 6) hA₀nn (by norm_num) (by norm_num) hγ1 hγ2
    have heq : (3 * |2 * γ| + 10) * (3 + 2 * (3 * A₀ + 4 * |γ| + |γ|))
        = (6 * |γ| + 10) * (3 + 6 * A₀ + 10 * |γ|) := by
      rw [abs_mul]; simp only [abs_two]; ring
    rw [heq]; exact h
  exact region_riemannZeta_R14_clean hγ hKbig hA₀nn (abs_nonneg γ) (abs_nonneg γ)
    (beta_lt_sigmaK hK0 hζ) hζ hBclean hA_pole hCclean

#print axioms riemannZeta_zero_free_region
