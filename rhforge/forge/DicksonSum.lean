import Mathlib.RingTheory.Polynomial.Chebyshev
import Mathlib.Tactic

/-!
# §1783 — Dickson identity: z^n + z^{-n} is a polynomial in z + z^{-1}

The atomic degree-reduction fact underlying the Chebyshev/Schur-Cohn on-circle route (with §1781/§1782).
The Dickson polynomials `D_n` (defined by `D₀=2, D₁=X, D_{n+2}=X·D_{n+1}−D_n`) satisfy, for every `z ≠ 0`,

    z^n + z^{-n}  =  D_n(z + z^{-1}).

Hence the *symmetric* Laurent monomial `z^n + z^{-n}` is a genuine polynomial in the single variable
`s = z + z^{-1}`.  This is the mechanism by which a self-inversive (palindromic) polynomial
`P(z) = Σ vᵢ zⁱ` of even degree collapses to a polynomial `Q(s)` of half the degree — turning "zeros of P
on the unit circle" into "roots of Q real in [-2,2]" via §1781/§1782, discharging §1777's owed input by
finite real-rootedness instead of the missing Herglotz theory.  EXIT0 / std-axioms.
-/

open Polynomial

variable {K : Type*} [Field K]

/-- Dickson polynomials (first kind, `a = 1`): `D₀ = 2`, `D₁ = X`, `D_{n+2} = X·D_{n+1} − D_n`. -/
noncomputable def dickson : ℕ → K[X]
  | 0 => C 2
  | 1 => X
  | (n + 2) => X * dickson (n + 1) - dickson n

/-- **Dickson sum identity.**  For `z ≠ 0`, `D_n(z + z⁻¹) = zⁿ + (z⁻¹)ⁿ`; i.e. the symmetric Laurent
monomial `zⁿ + z⁻ⁿ` is the value at `z + z⁻¹` of a polynomial of degree `n`. -/
theorem dickson_eval_sum (z : K) (hz : z ≠ 0) (n : ℕ) :
    (dickson n).eval (z + z⁻¹) = z ^ n + (z⁻¹) ^ n := by
  induction n using Nat.twoStepInduction with
  | zero => simp only [dickson, eval_C, pow_zero]; norm_num
  | one => simp [dickson]
  | more k ih0 ih1 =>
    rw [dickson, eval_sub, eval_mul, eval_X, ih0, ih1]
    simp only [inv_pow]
    field_simp
    ring
