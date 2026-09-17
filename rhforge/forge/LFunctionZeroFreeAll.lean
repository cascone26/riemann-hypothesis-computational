import LFunctionZeroFreeClean
import LFunctionZeroFreeCleanQuad

set_option maxHeartbeats 1000000

/-!
# §1568 ★★★★★★ — THE de la Vallée Poussin zero-free region for EVERY nontrivial Dirichlet character

THE UNIFIED CAPSTONE.  Case split on `χ² = 1`: quadratic → §1567, non-quadratic → §1552.  For every `N` and every
nontrivial Dirichlet character `χ` mod `N`:

  **∃ T₀ c₀ > 0, ∀ γ β, T₀ ≤ |γ| → LFunction χ(β+iγ) = 0 → β ≤ 1 - c₀/log|γ|.**

This completes the classical quantitative zero-free-region program (`ζ` §1119→§1464, non-quadratic L §1465→§1552,
quadratic L §1553→§1567) — the full de la Vallée Poussin theorem for the Riemann zeta function and ALL Dirichlet
L-functions, machine-checked end-to-end.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem LFunction_zero_free_region_clean_all {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) :
    ∃ T₀ c₀ : ℝ, 0 < c₀ ∧ ∀ γ β : ℝ, T₀ ≤ |γ| →
      LFunction χ ((β : ℂ) + I * (γ : ℝ)) = 0 → β ≤ 1 - c₀ / Real.log |γ| := by
  by_cases h : χ ^ 2 = 1
  · exact LFunction_zero_free_region_clean_quad hχ h
  · exact LFunction_zero_free_region_clean hχ h

#print axioms LFunction_zero_free_region_clean_all
