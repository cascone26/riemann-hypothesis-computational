# The Chebyshev On-Circle Criterion — a complete, kernel-verified Lean formalization

**Status:** COMPLETE (2026-08-29). 9 axiom-clean Lean bricks (§1781–1789), part of the CS-engine chain
(§1771–1789, 19 bricks total). All compile against Mathlib with axiom set ⊆ {propext, Classical.choice,
Quot.sound}. Location: Olympus `~/rhforge/forge/`, mirror `~/projects/riemann/rhforge/forge/`.

## What it proves (unconditional, self-contained)

For a **palindromic** polynomial `P(z) = Σ_{i≤2m} vᵢ zⁱ` (`vᵢ = v_{2m−i}`) and any nonzero root `ρ`:

> `ρ` lies **on the unit circle** `|z|=1`  ⟺  its Chebyshev image `s = ρ + ρ⁻¹` is a **real** root
> of the half-degree polynomial `Q` with `|s| ≤ 2`  (i.e. `s ∈ [−2,2]`).

This **reduces "zeros on the unit circle" to real-rootedness of a half-degree polynomial in [−2,2]** — a
finite, elementary path that avoids the (mathlib-absent) Herglotz / trigonometric-moment theory that the
usual Carathéodory–Fejér proof routes through.

## The 9 bricks

| Brick | Statement |
|-------|-----------|
| §1781 `JoukowskiOnCircle` | `z + z⁻¹` real with `|z+z⁻¹| ≤ 2` ⟺ `|z| = 1` (z≠0) |
| §1782 `JoukowskiSurj` | every `w ∈ [−2,2]` is `z+z⁻¹` for some `z` on the circle |
| §1783 `DicksonSum` | `Dₙ(z+z⁻¹) = zⁿ + z⁻ⁿ` (Dickson recurrence, two-step induction) |
| §1784 `DicksonComboSymmetric` | `zᵐ·Σ cⱼ Dⱼ(z+z⁻¹) = Σ cⱼ(z^{m+j}+z^{m−j})` |
| §1785 `RootOnCircleChebyshev` | given `P(ρ)=ρᵐQ(ρ+ρ⁻¹)`: root-on-circle ⟺ Q-root real in [−2,2] |
| §1786 `PowPairDickson` | `z^{m+k}+z^{m−k} = zᵐ·D_k(z+z⁻¹)` (k≤m) |
| §1787 `PalindromicPairSum` | `Σ_{i<2m+1} f i = f m + Σ_{k<m}(f(m+1+k)+f(m−1−k))` (Finset split) |
| §1788 `PalindromicDicksonFactor` | palindromic `P = zᵐ·Q(z+z⁻¹)` (the factorization) |
| §1789 `PalindromicRootOnCircle` | **HEADLINE**: composes §1788+§1781 end-to-end |

## Honest scope

- This is **durable formalization of the UNCONDITIONAL** on-circle reduction. Mathlib has no
  "self-inversive-polynomial-zeros-on-circle" result — this is a genuine, potentially PR-worthy contribution
  in its own right (cf. the earlier general-Bessel module).
- It is **NOT an RH proof**, and does not claim to be. It formalizes the *engine* of the Connes–van
  Suijlekom mechanism, not its RH-hard hypotheses.
- **The one remaining input** to close the full Carathéodory–Fejér engine: *Q has all real roots in [−2,2]*
  for the PSD-singular-Toeplitz kernel polynomial (the classical Schur–Cohn / Hermite–Bezoutian content).
  Mathlib lacks Sturm/Bezoutian/Hankel-signature machinery, so this is a large from-scratch build — and it is
  still unconditional (not where RH-hardness lives).
- **RH-hardness** lives in the Connes program's two open steps, both mapped RH-equivalent this session:
  (i) the Weil-form minimum eigenvalue being simple + even for all windows; (ii) the convergence θ_x → k_λ.

## Provenance
Full cycle-by-cycle derivation: `THE_BRIDGE.md` cycles 118–143 (route reframing, the plain-form artifact
retraction at cycle 129, the correct derivative-form object, and this Chebyshev sub-project 134–143).

## Update 2026-08-29 — wired to the real-rootedness module (§1812–1817)

The on-circle criterion now connects end-to-end to the Hankel-PSD iffs (§1806/§1811):
- **§1812** `PalindromicAllRoots` — all nonzero P-roots on `|z|=1` ⟺ all Chebyshev images `ρ+ρ⁻¹` real in `[−2,2]` (root-set lift of §1789).
- **§1813** `QpolyDef` — `Qpoly` as a genuine `Polynomial ℂ` (Dickson combination); `.roots` now available.
- **§1814** `PalindromicQpolyRoot` — `ρ≠0` is a P-root ⟺ `ρ+ρ⁻¹` is a `Qpoly`-root (Joukowski root correspondence).
- **§1815** `DicksonMapConj` — Dickson coefficients are conjugation-fixed (real).
- **§1816** `QpolyMapConj` — real `v` ⟹ `Qpoly` real-coefficiented ⟹ `Qpoly.roots` conjugation-closed.
- **§1817** `QpolyHankelIff` (CAPSTONE) — real `v`: complex Hankel form over `Qpoly.roots` PSD ⟺ all those roots real.

Chain: **P-roots-on-circle ↔ Qpoly-roots-real-in-[−2,2] ↔ Hankel-PSD**. §1817 realizes the core (real ↔ PSD);
the interval refinement (in-[−2,2] via §1811) and the full 2-to-1 `P.roots ↔ Qpoly.roots` multiset counting are
the remaining plumbing. This is the CS/Chebyshev mechanism's positivity condition, in Lean, axiom-clean.
