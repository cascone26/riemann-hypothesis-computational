import DLVPCharRegionScalar
import HpoleRiemannZeta

/-!
# §1470 — DLVP L-region with the ζ pole term discharged (`P → 1/(σ−1) + A₀`)

The ζ pole contribution `P = Re(−logDeriv ζ σ)` is the *same* object in both the ζ and the Dirichlet-L lanes, so
the proven ζ brick §1405 `hpole_riemannZeta` (`Re(−logDeriv ζ σ) ≤ 1/(σ−1) + ‖logDeriv (LFunctionTrivChar₁ 1) σ‖`)
discharges it here verbatim.  Feeding it into §1469 leaves only the two `L`-side bounds (`Q` at a zero of
`L(·,χ)`, `R` for `L(·,χ²)`) undischarged:

  `0 ≤ 3·(1/(σ−1) + A₀) + 4Q + R`,   `A₀ := ‖logDeriv (LFunctionTrivChar₁ 1) σ‖`.

This is the L-analogue of the ζ region assembly after its pole step.  EXIT0 / std-axioms.
-/

open Complex ArithmeticFunction LSeries DirichletCharacter
open scoped LSeries.notation

theorem dlvp_char_region_pole {N : ℕ} (χ : DirichletCharacter ℂ N) {σ t Q R : ℝ} (hσ : 1 < σ)
    (hzero : (-logDeriv (LSeries (↗χ)) ((σ : ℂ) + I * (t : ℝ))).re ≤ Q)
    (hchi2 : (-logDeriv (LSeries (↗(χ ^ 2))) ((σ : ℂ) + I * ((2 * t : ℝ) : ℂ))).re ≤ R) :
    0 ≤ 3 * (1 / (σ - 1) + ‖logDeriv (LFunctionTrivChar₁ 1) ((σ : ℂ))‖) + 4 * Q + R :=
  dlvp_char_region_scalar χ hσ (hpole_riemannZeta hσ) hzero hchi2

#print axioms dlvp_char_region_pole
