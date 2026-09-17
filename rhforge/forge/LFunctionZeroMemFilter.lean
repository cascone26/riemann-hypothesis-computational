import LFunctionZeroMemSupport

/-!
# §1496 — the DLVP zero lies in the **filtered** (half-disk) L support (`χ ≠ 1`; the region `hρ₀`)

L-analogue of ζ §1407.  The region composition §1336's `hρ₀` needs the kept zero `ρ = β+iγ` in the filtered
support `S = support.filter p` (the half-disk zeros).  From §1495 it is in the support, and if `p ρ` holds
(`ρ ∈ closedBall (σK+iγ) (R/2)`, i.e. `σK-β ≤ R/2`) then it is in the filter.  EXIT0 / std-axioms.
-/

open Metric Classical DirichletCharacter

theorem LFunction_zero_mem_filter {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {c : ℂ} {R' : ℝ} {ρ z₀ : ℂ} (hρ : ρ ∈ ball c R') (hLρ : LFunction χ ρ = 0)
    (hz₀ : z₀ ∈ ball c R') (hz₀ne : LFunction χ z₀ ≠ 0) (p : ℂ → Prop) (hp : p ρ) :
    ρ ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter p :=
  Finset.mem_filter.mpr ⟨LFunction_zero_mem_divisor_support hχ hρ hLρ hz₀ hz₀ne, hp⟩

#print axioms LFunction_zero_mem_filter
