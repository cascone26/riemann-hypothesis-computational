# Lean 4 Mathlib Audit for RH Formalization

**Date**: 2026-03-30
**Purpose**: Assess what Lean 4's Mathlib currently provides for Riemann Hypothesis-related formalization.

---

## 1. Complex Analysis in Mathlib

### Riemann Zeta Function -- YES, fully formalized

As of March 2025, the Riemann zeta function has been **fully formalized** in Mathlib, including analytic continuation and the functional equation. This was done by David Loeffler and published in the Annals of Formalized Mathematics (received 2025-03-06, accepted 2025-06-23).

Key definitions in `Mathlib.NumberTheory.LSeries.RiemannZeta`:

| Lean Name | Description |
|---|---|
| `riemannZeta` | The Riemann zeta function zeta : C -> C |
| `completedRiemannZeta` | Lambda(s) = pi^(-s/2) * Gamma(s/2) * zeta(s) |
| `completedRiemannZeta₀` | Entire function: Lambda_0(s) = Lambda(s) + 1/(s-1) - 1/s |
| `RiemannHypothesis` | Formal statement: all non-trivial zeros have Re(s) = 1/2 |

Key theorems proved:

- `riemannZeta_one_sub` -- Functional equation relating zeta(1-s) and zeta(s)
- `completedRiemannZeta_one_sub` -- Lambda(1-s) = Lambda(s)
- `riemannZeta_two` -- zeta(2) = pi^2/6 (Basel problem)
- `riemannZeta_eulerProduct` -- Euler product formula
- `riemannZeta_ne_zero_of_one_le_re` -- zeta(s) != 0 for Re(s) >= 1
- `differentiableAt_riemannZeta` -- Differentiable everywhere except s = 1
- `riemannZeta_residue_one` -- Residue at s=1 equals 1
- `riemannZeta_neg_two_mul_nat_add_one` -- Trivial zeros at negative even integers
- `zeta_eq_tsum_one_div_nat_cpow` -- Dirichlet series representation for Re(s) > 1

**Approach**: Uses theta functions and Mellin transforms (not contour integrals). The construction proceeds:
1. Jacobi theta function theta(tau) = sum of e^(pi*i*n^2*tau)
2. Poisson summation to get the transformation law
3. Mellin transform connects theta to completed zeta function
4. Analytic continuation follows from the integral representation

Source: https://arxiv.org/abs/2503.00959
Docs: https://leanprover-community.github.io/mathlib4_docs/Mathlib/NumberTheory/LSeries/RiemannZeta.html

### Gamma Function -- YES

Fully formalized in `Mathlib.Analysis.SpecialFunctions.Gamma.*` (Basic, Beta, BohrMollerup, Deligne, Deriv, Digamma).

| Lean Name | Description |
|---|---|
| `Complex.Gamma` | Gamma function for complex variable |
| `Real.Gamma` | Gamma function for real variable |

Properties proved:
- Euler integral representation for Re(s) > 0
- Recurrence: Gamma(s+1) = s * Gamma(s)
- Gamma(n+1) = n! for naturals
- Meromorphic continuation to all of C (via recurrence)
- Conjugate symmetry
- Positivity for positive reals
- Convention: Gamma(0) = 0, Gamma(-n) = 0 (junk values at poles)

Docs: https://leanprover-community.github.io/mathlib4_docs/Mathlib/Analysis/SpecialFunctions/Gamma/Basic.html

### Analytic Continuation -- YES (for specific functions)

Mathlib does NOT have a general "analytic continuation" procedure as an abstract operation. Instead, specific functions (zeta, Dirichlet L-functions, Gamma) are constructed with their continuations built in.

The principle of analytic continuation (identity theorem) IS formalized:
- If two holomorphic functions agree on a set with an accumulation point, they agree everywhere on the connected domain.

### Dirichlet Series / L-functions -- YES

Formalized in `Mathlib.NumberTheory.LSeries.*`:

- `Mathlib.NumberTheory.LSeries.DirichletContinuation` -- Analytic continuation of Dirichlet L-functions
- `DirichletCharacter.LFunction` -- Unique meromorphic function agreeing with the L-series
- `DirichletCharacter.completedLFunction` -- Completed L-function with gamma factors
- Functional equations for primitive characters
- Non-vanishing on Re(s) >= 1
- **Dirichlet's theorem on primes in arithmetic progressions** (proved!)

Docs: https://leanprover-community.github.io/mathlib4_docs/Mathlib/NumberTheory/LSeries/DirichletContinuation.html

### Meromorphic Functions -- YES (recently added)

Mathlib now has `Mathlib.Analysis.Meromorphic.*` modules:

| Module | Content |
|---|---|
| `Meromorphic.Basic` | `MeromorphicAt`, `MeromorphicOn` predicates |
| `Meromorphic.NormalForm` | `MeromorphicNFAt`, `MeromorphicNFOn` |
| `Meromorphic.Divisor` | Divisor of a meromorphic function (order at each point) |
| `Meromorphic.IsolatedZeros` | Principle of isolated zeros for meromorphic functions |
| `Meromorphic.Complex` | Complex-specific meromorphic function theory |

Note: The "undergrad_todo" page and "mathlib-overview" page still list meromorphic functions as missing, but the actual codebase now has them. These overview pages appear to be outdated.

Docs: https://leanprover-community.github.io/mathlib4_docs/Mathlib/Analysis/Meromorphic/NormalForm.html

### Complex Integration / Contour Integrals -- PARTIAL

**Present**:
- Cauchy integral formula (circle integrals)
- Liouville's theorem
- Maximum modulus principle
- Schwarz lemma
- Removable singularity theorem

**Missing** (still listed on undergrad_todo):
- General contour integrals along paths
- Winding number
- Laurent series
- Residue theorem

The zeta function formalization deliberately avoids contour integrals by using the theta function / Mellin transform approach instead.

---

## 2. Arithmetic Functions in Mathlib

### Divisor Function sigma(n) -- YES

Defined in `Mathlib.NumberTheory.ArithmeticFunction.*`:

```
ArithmeticFunction.sigma k n = sum over divisors d of n: d^k
```

For k=1, this gives the standard sum-of-divisors function sigma_1(n).

The sigma function is proved to be multiplicative. Various identities involving divisor sums are available.

Docs: https://leanprover-community.github.io/mathlib4_docs/Mathlib/NumberTheory/ArithmeticFunction.html

### Euler's Totient Function -- YES

Defined as `Nat.totient n` in `Mathlib.Data.Nat.Totient`:

```
n.totient = card of {a in range(n) | coprime(n, a)}
```

Properties proved:
- Divisor sum formula: n = sum of phi(d) over divisors d of n
- `totient_mul` for coprime arguments
- `totient_prime_pow` for prime powers

Docs: https://leanprover-community.github.io/mathlib4_docs/Mathlib/Data/Nat/Totient.html

### Mobius Function -- YES

Defined as `ArithmeticFunction.moebius` (notation: mu) in `Mathlib.NumberTheory.ArithmeticFunction.Moebius`:

- mu(n) = (-1)^k if n is squarefree with k distinct prime factors
- mu(n) = 0 if n is not squarefree

Properties proved:
- `isMultiplicative_moebius` -- multiplicativity
- `moebius_apply_prime` -- mu(p) = -1
- `moebius_apply_prime_pow` -- mu(p^k) behavior
- **Mobius inversion** -- `sum_eq_iff_sum_mul_moebius_eq`

### Euler-Mascheroni Constant -- YES

Defined as `Real.eulerMascheroniConstant` in `Mathlib.NumberTheory.Harmonic.EulerMascheroni`:

Properties proved:
- `Real.tendsto_harmonic_sub_log` -- harmonic(n) - log(n) converges to gamma
- `Real.one_half_lt_eulerMascheroniConstant` -- gamma > 1/2
- `Real.eulerMascheroniConstant_lt_two_thirds` -- gamma < 2/3
- Strict monotonicity of approximating sequences
- Bounds: 1/2 < gamma < 2/3

Docs: https://leanprover-community.github.io/mathlib4_docs/Mathlib/NumberTheory/Harmonic/EulerMascheroni.html

### Harmonic Numbers -- YES

Defined as `harmonic n` in `Mathlib.NumberTheory.Harmonic.Defs`:

```
harmonic n = sum over i in range(n): 1/(i+1)
```

Returns a rational number. Basic recurrence (`harmonic_succ`) and base case (`harmonic_zero`) proved.

---

## 3. Existing RH-Related Formalization Work

### Statement of RH in Lean 4 -- YES, in Mathlib itself

`RiemannHypothesis` is defined in Mathlib. The statement says: all zeros of zeta(s), other than the trivial zeros at negative even integers, have Re(s) = 1/2.

This is the **official** formalization, now part of the Mathlib library.

### Kontorovich/Gomes Lean-RH -- YES (older, standalone)

Repository: https://github.com/AlexKontorovich/Lean-RH

Uses the Dirichlet eta function formulation: if 0 < Re(s) < 1 and DirichletEta(s) = 0, then Re(s) = 1/2. Mostly standalone (not depending on Mathlib), with a separate `impl` file showing the axioms can be derived from Mathlib.

### Millennium Prize Problems -- YES

Repository: https://github.com/lean-dojo/LeanMillenniumPrizeProblems

Formalizes the Clay Mathematics Institute statement of RH as `Millennium.RiemannHypothesis`, with an equivalence lemma to Mathlib's `RiemannHypothesis`.

### Robin's Inequality -- NO formalization exists

No Lean formalization of Robin's inequality (sigma(n) < e^gamma * n * ln(ln(n)) for n > 5040) has been found.

### Li's Criterion -- NO formalization exists

No Lean formalization of Li's criterion has been found.

### Other Equivalent Formulations -- NO

No other equivalent formulations of RH (Lagarias inequality, Redheffer matrix, etc.) appear to be formalized in Lean.

---

## 4. Feasibility Assessment

### Option A: Robin's Inequality Statement

**sigma(n) < e^gamma * n * ln(ln(n)) for all n > 5040, assuming RH**

What exists:
- [x] ArithmeticFunction.sigma (divisor sum function)
- [x] Real.eulerMascheroniConstant (gamma, with bounds 1/2 < gamma < 2/3)
- [x] Real.exp (exponential function)
- [x] Real.log (natural logarithm)
- [x] Natural number arithmetic, inequalities
- [x] RiemannHypothesis (the statement to assume or relate to)

What needs to be built:
- [ ] Statement connecting RH to Robin's inequality (the actual theorem)
- [ ] Possibly better numerical bounds on gamma (current: 1/2 < gamma < 2/3; may need tighter for computational verification of small cases)

**Verdict: MOST FEASIBLE.** All the building blocks exist in Mathlib. The statement itself only requires combining existing definitions with basic arithmetic. You would write something like:

```lean
theorem robin_inequality_iff_RH :
  RiemannHypothesis ↔ ∀ n : ℕ, n > 5040 →
    (ArithmeticFunction.sigma 1 n : ℝ) <
      Real.exp Real.eulerMascheroniConstant * n * Real.log (Real.log n) := sorry
```

This is a clean statement that uses only existing Mathlib definitions. PROVING it is another matter entirely (that is a deep theorem by Robin, 1984), but STATING it is immediately achievable.

### Option B: Statement of RH

**Already done.** `RiemannHypothesis` exists in Mathlib. No work needed for the statement. If the goal is to state an equivalent formulation, that requires proving equivalence.

### Option C: Li's Criterion

**The Li coefficients lambda_n = sum over non-trivial zeros rho of (1 - (1 - 1/rho)^n)**

What's missing:
- No formalization of "the set of non-trivial zeros of zeta"
- No way to sum over zeros (would need to enumerate them or use a suitable limit)
- The concept of "zeros of a meromorphic function" exists via the divisor/order machinery, but extracting them as a countable set to sum over is nontrivial

**Verdict: HARD.** Would require substantial new infrastructure.

### Recommended Path

1. **Immediate win**: State Robin's inequality as a Lean 4 theorem (the statement, not the proof). All ingredients exist.
2. **Medium-term**: Prove Robin's inequality for specific small cases (n <= some bound) computationally.
3. **Long-term**: The equivalence proof (Robin's inequality <-> RH) is a serious research-level formalization project.

---

## 5. Lean 4 + Mathlib Installation on macOS (Apple Silicon)

### Prerequisites

- Git and curl (should already be installed)
- Homebrew (recommended)

### Step 1: Install elan (Lean version manager)

```bash
curl https://elan.lean-lang.org/elan-init.sh -sSf | sh
source ~/.elan/env
```

elan manages Lean versions per-project (like rustup for Rust or nvm for Node).

### Step 2: Install VS Code + Lean 4 extension

```bash
brew install --cask visual-studio-code
```

Then install the `leanprover.lean4` extension from the VS Code marketplace.

### Step 3: Create a new project with Mathlib

```bash
lake new riemann-lean math
```

The `math` keyword at the end adds mathlib4 as a dependency automatically. This creates:
- `lakefile.lean` -- build configuration with Mathlib dependency
- `lean-toolchain` -- specifies the Lean version
- Git repository initialized

### Step 4: Build / fetch Mathlib cache

```bash
cd riemann-lean
lake exe cache get    # downloads pre-built Mathlib oleans (MUCH faster than building)
lake build
```

Building Mathlib from source takes hours. The cache download takes minutes.

### Step 5: Verify

Create a file and add:
```lean
import Mathlib

#check RiemannHypothesis
#check ArithmeticFunction.sigma
#check Real.eulerMascheroniConstant
```

### Notes

- Apple Silicon is fully supported. Lean 4 has native ARM64 builds.
- Mathlib is large (~4GB with cache). Ensure sufficient disk space.
- The `lean-toolchain` file pins the exact Lean version. Don't change it independently of Mathlib.
- `lake update` + `lake exe cache get` to update Mathlib.

Official guide: https://leanprover-community.github.io/install/macos.html (redirects to https://lean-lang.org/install/manual/)
Project setup: https://leanprover-community.github.io/install/project.html

---

## Summary Table

| Component | In Mathlib? | Lean Name | Notes |
|---|---|---|---|
| Riemann zeta function | YES | `riemannZeta` | Full analytic continuation + functional equation |
| Gamma function | YES | `Complex.Gamma` / `Real.Gamma` | Meromorphic, recurrence, factorial |
| Analytic continuation | YES (specific) | Built into zeta/L-function defs | No general abstract procedure |
| Dirichlet L-functions | YES | `DirichletCharacter.LFunction` | Continuation + functional equation |
| Meromorphic functions | YES | `MeromorphicAt` / `MeromorphicOn` | Recently added |
| Contour integrals | NO | -- | Cauchy integral formula exists, but not general contours |
| Residue theorem | NO | -- | Listed as missing in undergrad_todo |
| Divisor function sigma | YES | `ArithmeticFunction.sigma` | sigma k n = sum of d^k over divisors |
| Euler totient | YES | `Nat.totient` | With divisor sum formula |
| Mobius function | YES | `ArithmeticFunction.moebius` | With Mobius inversion |
| Euler-Mascheroni gamma | YES | `Real.eulerMascheroniConstant` | Bounds: 1/2 < gamma < 2/3 |
| Harmonic numbers | YES | `harmonic` | Rational-valued, basic properties |
| RH statement | YES | `RiemannHypothesis` | In Mathlib, zeros have Re = 1/2 |
| Robin's inequality | NO | -- | Not formalized anywhere |
| Li's criterion | NO | -- | Would need zero-enumeration infrastructure |

---

## Key Sources

- Paper: "Formalizing zeta and L-functions in Lean" (Loeffler, 2025) -- https://arxiv.org/abs/2503.00959
- Mathlib RiemannZeta docs -- https://leanprover-community.github.io/mathlib4_docs/Mathlib/NumberTheory/LSeries/RiemannZeta.html
- Lean-RH (Kontorovich/Gomes) -- https://github.com/AlexKontorovich/Lean-RH
- LeanMillenniumPrizeProblems -- https://github.com/lean-dojo/LeanMillenniumPrizeProblems
- Mathlib overview -- https://leanprover-community.github.io/mathlib-overview.html
- Mathlib undergrad todo -- https://leanprover-community.github.io/undergrad_todo.html
