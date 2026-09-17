# The Hermite Real-Rootedness Module — a from-scratch Lean formalization

**Status:** COMPLETE (2026-08-29) — Hermite criterion as a single iff `hankel_psd_iff_all_real` (conj-closed, all-simple roots). 28 files, all axiom-clean — real-rootedness iff (§1806) + interval iff (§1811) + on-circle wiring §1812–1817 (Qpoly-level Hankel criterion connecting on-circle ⊗ real-rootedness). 14 module files / 15 theorems (§1790–1803), all compiling against Mathlib
with axiom set ⊆ {propext, Classical.choice, Quot.sound}. Location: Olympus `~/rhforge/forge/`, mirror
`~/projects/riemann/rhforge/forge/`.

## Why this module

The on-circle criterion (§1781–1789, `ON_CIRCLE_CRITERION.md`) reduces "palindromic-polynomial zeros on the
unit circle" to "the half-degree polynomial `Q` has all real roots in `[−2,2]`". Closing the Connes–van
Suijlekom / Carathéodory–Fejér engine therefore needs a **real-rootedness criterion**: a decidable positivity
condition equivalent to "all roots real (and in an interval)". The classical answer is **Hermite's theorem**:
for a real polynomial, the number of *distinct real roots* equals the **signature** of the Hankel matrix of
Newton power sums `s_{i+j} = ∑_roots r^{i+j}`. Mathlib has **none** of this (verified: no Sturm, no
Bezoutian/Hankel-signature; `Hermite`=Hermite polynomials, `Bezout`=domain algebra). This is a genuine,
PR-worthy, from-scratch build — and, crucially, it is **unconditional**: RH-hardness does NOT live here (it
lives in the Weil ∀-scales positivity, cf. `THE_BRIDGE.md` cycles 154–155).

## The 11 bricks

| Brick | Module | Statement |
|-------|--------|-----------|
| §1790 | `HankelPowerSumSOS` | family `r:Fin n→ℝ`: `∑_{i,j}(∑_k r_k^{i+j})x_i x_j = ∑_k(∑_i x_i r_k^i)²` (≥0) — **all-real ⟹ Hankel PSD** |
| §1791 | `WeightedHankelInterval` | weighted `w`: same with weight; PSD when `w≥0`; **weight `4−r²≥0` ⟹ [−2,2] interval-localization** |
| §1792 | `MultisetHankel` | lift to `R:Multiset ℝ` (single-root lemma + Multiset induction) — bridge to `Polynomial.roots` |
| §1793 | `PolynomialHankel` | Polynomial-level: root-power-sum Hankel form of any `p:ℝ[X]` is PSD (Hermite necessary condition); `s₀=#roots` |
| §1794 | `ConjPairNeg` | **converse SEED**: non-real `r` ⟹ real coeffs make conj-pair contribute `i²+(-i)²=−2 <0` |
| §1795 | `ComplexHankel` | complex-root form: `∑_{i,j}(∑_{r∈R}r^{i+j})c_i c_j = ∑_{r∈R}(∑_i c_i r^i)²` over ℂ, real `c` |
| §1796 | `ComplexHankelReal` | conj-closed `R` ⟹ complex Hankel form is **real** (im=0) — makes "PSD" a meaningful hypothesis |
| §1797 | `HankelConverseCore` | **converse CORE**: a real vector isolating a non-real root drives the form to `-2` (not PSD) |
| §1798 | `RealLinearSurj` | at a non-real node a real linear form hits *every* complex value (interpolation engine, generalizes §1794) |
| §1799 | `HankelPairConverse` | full converse, single-pair case (§1798∘§1797): a non-real conjugate pair's Hankel form is not PSD |
| §1800 | `ConjClosedProdReal` | conj-closed `S` ⟹ `∏(X−C s)` has real coefficients (isolating product is real) |
| §1801 | `ProdEvalNeZero` | `∏_{s∈S}(r₀−s) ≠ 0` for `r₀∉S` (isolating product nonzero at the target) |
| §1802 | `RealPolyCoeffVector` | real-coeff `P:ℂ[X]` ⟹ real vector `c` with `∑ c_i r^i = P.eval r` (coeff bridge) |
| §1803 | `HankelFullConverse` | full Hermite converse: conj-closed roots with a non-real element ⟹ Hankel form not PSD (`= -2`) |
| §1804 | `ComplexHankelNonnegReal` | all roots real ⟹ complex Hankel form PSD (the ⟸ direction) |
| §1805 | `ConjClosedDecomp` | conj-closed R + simple non-real r₀ ⟹ R = r₀::conj r₀::rest, rest conj-closed, r₀,conj r₀∉rest |
| §1806 | `HankelIff` | **CAPSTONE — single iff**: conj-closed all-simple R, complex Hankel form PSD ⟺ all roots real |

## Logical structure

- **Necessary direction (COMPLETE):** all roots real ⟹ Hankel PSD (§1790, §1792, §1793); roots in [−2,2] ⟹
  secondary `(4−r²)`-Hankel PSD (§1791). These give the "easy" half at family / weighted / multiset /
  polynomial levels.
- **Converse (hard direction, IN PROGRESS):** PSD ⟹ all roots real. Pieces in place: the **core logic** (§1797: an isolating real vector forces the form to `-2`), the **object** it acts
  on (§1795, complex Hankel form = `∑_r P_c(r)²`), its **reality** for conj-closed root sets (§1796), and the
  **seed** (§1794, a non-real root's conjugate pair can be driven to `−2`).
  - **Remaining (the one real gap):** the **isolation** — for a non-real root `r₀∈R`, construct a *real*
    polynomial `P` with `P(r₀)=i` and `P=0` on all other roots. Then (real coeffs ⟹ `P(r̄₀)=−i`
    automatically) the complex Hankel form evaluates to `i²+(−i)²+0 = −2 < 0`, contradicting PSD ⟹ all roots
    real. Construction: `P = Q·L` with `Q = ∏_{other roots}(t−s)` (real, conj-closed product) and `L` a real
    linear map hitting `L(r₀)=i/Q(r₀)` (possible since `r₀` non-real, cf. §1794's own construction). Mathlib
    has `Lagrange.interpolate`. **§1799 closes this for the single-pair case** (empty "rest") end-to-end;
    the general case adds only the real product `Q = ∏_{other roots}(t−s)` vanishing on the remaining roots
    (§1798 supplies the linear factor `L` with `L(r₀)=i/Q(r₀)`). Bounded from-scratch build, RH-hardness-free.

## Honest scope

- Genuinely novel formalization: Mathlib lacks any real-rootedness / Hankel-signature theory. Even the
  necessary direction (§1790–1793) is new.
- **NOT an RH proof.** This is the unconditional *engine* of the CS/on-circle mechanism, not its RH-hard
  hypotheses (the Weil ∀-scales positivity converse — field-wide open, cf. arXiv:2607.24830, `THE_BRIDGE.md`
  cycle 154).
- Independent July-2026 literature (arXiv:2607.24830, 30-digit numerical) confirms our analytic findings
  (knife-edge λ_a, μ₁ open constant, prime↔zero duality) and hits the same wall — but has **no**
  formalization; this module is the axis where we are ahead of the published frontier.

## Provenance
Full cycle-by-cycle derivation: `THE_BRIDGE.md` cycles 156–161. Necessary-direction foundation cycle 156;
interval cycle 157; multiset/polynomial bridge cycles 158–159; converse seed/object/reality cycles 159–161.
