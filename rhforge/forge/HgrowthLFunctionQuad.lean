import LogDerivLFunctionTrivCharSplit
import ProdCorrectionReLe
import HgrowthRiemannZetaScalar

set_option maxHeartbeats 1000000

/-!
# §1556 — QUADRATIC (`χ²=1`) growth bound: `Re(-logDeriv (LFunction χ²) c) ≤ B + ∑_{p|N} log p`

Combines the three quadratic-case bricks.  When `χ²=1`, `LFunction (χ²) = LFunction 1 = LFunctionTrivChar N`, whose
logDeriv splits (§1553) as `logDeriv(∏_{p|N}(1-p^{-·})) + logDeriv ζ`.  Taking `Re(-·)` and summing:
`Re(-logDeriv(LFunction χ²) c) = Re(-logDeriv(∏) c) + Re(-logDeriv ζ c) ≤ (∑_{p|N} log p) + B`, where the
product term is bounded by §1554 and the ζ term by §1555 (the pole-aware BC scalar, `hball` at `c`).  This is the
replacement for §1540 in the quadratic case — the χ side of the region is unchanged (χ itself is nontrivial).
EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric

theorem hgrowth_LFunction_quad {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ2eq : χ ^ 2 = 1)
    {c : ℂ} {R R' Ms B : ℝ} (hR : 0 < R) (hRR' : R < R') (hcre : 1 < c.re)
    (hball : (1 : ℂ) ∉ ball c R')
    (hζsphere : ∀ z ∈ sphere c R, ‖riemannZeta z‖ ≤ Ms)
    (hByscalar : 4 * Real.log ((Ms + 1) * (3 + 1 / (c.re - 1))) / (R / 2) ≤ B) :
    (-logDeriv (LFunction (χ ^ 2)) c).re ≤ B + ∑ p ∈ N.primeFactors, Real.log p := by
  have hLeq : LFunction (χ ^ 2) = LFunctionTrivChar N := by rw [hχ2eq]
  rw [hLeq]
  have hsplit := logDeriv_LFunctionTrivChar_split (N := N) hcre
  have hre : (-logDeriv (LFunctionTrivChar N) c).re
      = (-logDeriv (fun w => ∏ p ∈ N.primeFactors, (1 - (p : ℂ) ^ (-w))) c).re
        + (-logDeriv riemannZeta c).re := by
    rw [hsplit, neg_add, Complex.add_re]
  rw [hre]
  have hcorr := prod_correction_re_le (N := N) hcre
  have hζg := hgrowth_riemannZeta_scalar hR hRR' hcre hball hζsphere hByscalar
  linarith [hcorr, hζg]

#print axioms hgrowth_LFunction_quad
