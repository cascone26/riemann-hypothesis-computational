# RiemannLean — verified formalization manifest

Machine-checked Lean 4 formalization of the Riemann sign theorem's skeleton and the start of the
Weil positivity criterion. **All theorems below build clean (`lake build RiemannLean`) and depend
only on the trusted kernel axioms `{propext, Classical.choice, Quot.sound}` — 0 sorries.**
Verified via the propose→verify loop (`verify_loop.py`), whose gate enforces the axiom set via
`#print axioms` (LEDGER §733/§735). Provenance: `handles.jsonl`, LEDGER §743–751.

## Modules

### `RiemannLean/Bessel.lean` — Bessel `J₀`, `J₁` and positivity (mathlib has no Bessel functions)
- `besselJ0`, `besselJ1` — definitions via power series
- `besselJ0_series_summable`, `besselJ0_zero` (`J₀(0)=1`)
- **`besselJ0_pos`** — `J₀(x) > 0` for `|x| < 2` (alternating-series bound `J₀ ≥ 1−(x/2)²`)
- **`besselJ1_pos`** — `J₁(x) > 0` for `0<x`, `x²<8`
- `besselJ0_pos_inv_sqrt`, `besselJ1_pos_inv_sqrt` — `> 0` at `1/√n`, `n ≥ 2`
- `besselJ0_pos_of_prime`, `besselJ1_pos_of_prime` — stated for primes
- `besselJ0_prod_pos` — `∏ J₀(1/√p) > 0`

### `RiemannLean/SignTheorem.lean` — the sign theorem's logical skeleton
- `onset_sign_eq_moebius` — for squarefree `n`, `P>0`: `sgn(μ(n)·P) = μ(n)`
- `bessel_ratio_prod_pos` — `∏ J₁/J₀(1/√p) > 0`
- **`sign_theorem_onset`** — for squarefree `n`, IF the onset `D = μ(n)·([∏J₁/J₀]·[∏J₀])`
  (hypothesis `hD` = the §694 Jacobi–Anger onset factorization; its resummation §714/§720 is the
  one open analytic step), THEN `sgn D = μ(n)`.

The elementary Bessel positivity (the sign theorem's core, asserted for months) is now proven;
the single unproven analytic input is an **explicit labeled hypothesis**, not hidden.

### `RiemannLean/WeilForm.lean` — start of the RH positivity criterion
- `archWeight` — the archimedean density `Re[½(ψ(¼+ir/2)+ψ(¼−ir/2))]−log π` (on `Complex.digamma`)
- `archWeight_even` — even in `r`

### `RiemannLean/TailLemma.lean`, `TailAnalytic.lean` — Route-6 tail kernel (pre-existing, 0 sorry)
- `tail_kernel`, `tail_F_nonneg` — `cos x · cosh w ≤ 1` under the tail hypotheses.

## Honest scope
No Handle for RH. This verifies **consequence-side** ingredients (the sign theorem's elementary
core) and **begins** the criterion side (Weil archimedean weight). RH itself is the positivity of
the full Weil form, whose hard core is the arithmetic/prime fluctuation (LEDGER §749) — the open
frontier. What is genuinely novel here: no prior formalization of the sign theorem or of Bessel
`J₀/J₁` positivity existed; this exploits the automated verified-formalization loop.

## Next verified increments (loaded)
- Bessel `J₀` monotone decreasing on `[0, first zero]` (`J₀' = −J₁`) → uniform lower bound.
- Weil `poleTerm`, `primeTerm` (on `ArithmeticFunction.vonMangoldt`) → the finite Weil form.
- The finite Weil positivity certificate structure (organ-2 companion).
