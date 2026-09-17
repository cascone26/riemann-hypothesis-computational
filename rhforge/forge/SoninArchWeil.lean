import Mathlib

open Matrix

namespace SoninArchWeil

variable {n : ℕ}

/-- **Connes 2026 §1195 skeleton — archimedean Weil positivity from the Sonin lower bound.**

`S` is the Sonin-space orthogonal projection `𝔖` (a Hermitian idempotent); `V` is the scaling
operator `ϑ(g)`.  Connes–Consani's semilocal-trace-formula lower bound

  `Tr(V S Vᴴ) ≤ Winf`   (with `Winf = W_∞(g*g*)`)

is their *actively-researched analytic input* — Report-tier, taken here as an explicit
hypothesis, NOT proven inside this brick.  Given that bound, archimedean Weil positivity
`0 ≤ Winf` follows.

The content that is machine-verified here (the HANDLE) is the positivity of the Sonin trace
term `0 ≤ Tr(V S Vᴴ)`:
  * an orthogonal projection `S` is positive-semidefinite (`S = S·Sᴴ`),
  * PSD is preserved under conjugation `A ↦ V A Vᴴ`,
  * a PSD matrix has nonnegative trace.
This is exactly the reason C-C's bound implies positivity; the bound itself remains C-C's
open input. -/
theorem arch_weil_nonneg_of_sonin_bound
    (V S : Matrix (Fin n) (Fin n) ℝ)
    (hsym : Sᴴ = S) (hidem : S * S = S)
    (Winf : ℝ) (hbound : (V * S * Vᴴ).trace ≤ Winf) :
    0 ≤ Winf := by
  have hS : S.PosSemidef := by
    have h := Matrix.posSemidef_self_mul_conjTranspose S
    rwa [hsym, hidem] at h
  have hVSV : (V * S * Vᴴ).PosSemidef := hS.mul_mul_conjTranspose_same V
  have htr : 0 ≤ (V * S * Vᴴ).trace := hVSV.trace_nonneg
  linarith

#print axioms arch_weil_nonneg_of_sonin_bound

end SoninArchWeil
