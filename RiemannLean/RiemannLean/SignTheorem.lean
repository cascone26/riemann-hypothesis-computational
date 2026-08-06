import RiemannLean.Bessel
import Mathlib.NumberTheory.ArithmeticFunction.Moebius
import Mathlib.Data.Real.Sign

/-!
# The Sign Theorem — formally verified logical skeleton

The Riemann **Sign Theorem** (LEDGER §690–698, `SIGN_THEOREM.md`): for every squarefree `n`,
the θ-twisted alternating zero sum reads the Möbius function,
`sgn Re[e^{-ic} G(log n)] = μ(n)`.

Its mechanism (§693–694, §722): in the Selberg-oscillator / Jacobi–Anger picture the coefficient
of `n^{iγ}` onsets at order `ω(n)` as
`D(n) = (-1)^{ω(n)} · ∏_{p∣n} [J₁/J₀](1/√p) · ∏_q J₀(1/√q)`,
and `sgn D(n) = μ(n)` because **every Bessel factor is strictly positive** (`1/√p ≤ 1/√2` lies
below the first zeros of `J₀, J₁`) and `(-1)^{ω(n)} = μ(n)` for squarefree `n`.

This file verifies that mechanism **completely, except** for the one genuinely-hard analytic
input — that the onset *equals* `μ(n)` times the positive Bessel product (the Jacobi–Anger
resummation on the zero set, §714/§720, still open) — which appears here as an explicit
hypothesis `hD`.  Everything else is machine-checked (trusted kernel axioms only):

* the Bessel positivity `∏ J₀(1/√p) > 0` — `Bessel.besselJ0_prod_pos` (fully proven);
* the sign algebra `sgn(μ(n)·P) = μ(n)` for squarefree `n`, `P > 0` — `onset_sign_eq_moebius`.
-/

open ArithmeticFunction

/-- **Sign-theorem logical core.**  For squarefree `n`, `P > 0`: `sgn(μ(n)·P) = μ(n)`. -/
theorem onset_sign_eq_moebius {n : ℕ} (hn : Squarefree n) {P : ℝ} (hP : 0 < P) :
    Real.sign ((moebius n : ℝ) * P) = (moebius n : ℝ) := by
  have hne : moebius n ≠ 0 := moebius_ne_zero_iff_squarefree.mpr hn
  rcases moebius_ne_zero_iff_eq_or.mp hne with h1 | h1
  · have : (moebius n : ℝ) = 1 := by rw [h1]; norm_num
    rw [this, one_mul, Real.sign_of_pos hP]
  · have : (moebius n : ℝ) = -1 := by rw [h1]; norm_num
    rw [this, neg_one_mul, Real.sign_of_neg (by linarith)]

/-- **Sign Theorem (conditional on the onset factorization).**  For squarefree `n`, if the
diagonal onset `D` factors as `μ(n)` times a `J₀`-product over primes `p ≥ 2` — the analytic
input `hD` (Jacobi–Anger onset at order `ω(n)`, §722; its resummation is the one open step) —
then the onset's sign is exactly `μ(n)`.  The positivity of the Bessel product and the sign
algebra are fully proven; only `hD` is assumed. -/
theorem sign_theorem_conditional {n : ℕ} (hn : Squarefree n)
    {S : Finset ℕ} (hS : ∀ p ∈ S, 2 ≤ p)
    {D : ℝ} (hD : D = (moebius n : ℝ) * ∏ p ∈ S, besselJ0 (1 / Real.sqrt p)) :
    Real.sign D = (moebius n : ℝ) := by
  rw [hD]
  exact onset_sign_eq_moebius hn (besselJ0_prod_pos hS)
