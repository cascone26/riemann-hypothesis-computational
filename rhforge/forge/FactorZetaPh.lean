import Mathlib.Analysis.SpecialFunctions.Complex.Log

/-!
# §1328 — Trivial local factorization `ζ =ᶠ P·(ζ/P)` where `P ≠ 0`

Wherever `P` is nonvanishing on a neighbourhood of `z`, `ζ = P·(ζ/P)` there (`mul_div_cancel`).  With `h := ζ/P`
this supplies the factorization hypothesis `ζ =ᶠ P·h` that §1327 needs to decompose `logDeriv ζ`.  At the DLVP test
point `z = σ+iγ` (`Re > 1`, so `ζ(z) ≠ 0` and, `z` being no zero, `P(z) = ∏(z-ρ) ≠ 0`), `h = ζ/P` is differentiable
and nonvanishing at `z` — so §1327+§1328 give the decomposition `logDeriv ζ z = ∑1/(z-ρ) + logDeriv(ζ/P) z` cleanly.
(The disk-wide nonvanishing of `ζ/P` — for the §1325 Borel–Carathéodory bound on `logDeriv(ζ/P)` — is the separate
multiplicity/divisor construction.)  EXIT0 / std-axioms.
-/

open Complex Filter Topology

theorem eventuallyEq_mul_div_of_ne {ζ P : ℂ → ℂ} {z : ℂ} (hP : ∀ᶠ s in 𝓝 z, P s ≠ 0) :
    ζ =ᶠ[𝓝 z] fun s => P s * (ζ s / P s) := by
  filter_upwards [hP] with s hs
  field_simp

#print axioms eventuallyEq_mul_div_of_ne
