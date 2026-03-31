# Hilbert-Polya Operator Candidates: Comprehensive Survey

## Background

The Hilbert-Polya conjecture (attributed informally to Hilbert and Polya, ~1910s-1930s) states that the nontrivial zeros of the Riemann zeta function correspond to eigenvalues of a self-adjoint operator. If such an operator exists, the Riemann Hypothesis follows immediately from self-adjointness (which guarantees real eigenvalues). The conjecture was dramatically strengthened by Montgomery's pair correlation conjecture (1973) and subsequent numerical work by Odlyzko, showing that zeta zero statistics match the Gaussian Unitary Ensemble (GUE) of random matrix theory — implying the hypothetical operator should resemble a quantum Hamiltonian of a chaotic system without time-reversal symmetry.

This document surveys every major proposed operator candidate in the literature.

---

## 1. Berry-Keating Conjecture (1999)

### Key References
- Berry & Keating, "H = xp and the Riemann Zeros," in *Supersymmetry and Trace Formulae* (Kluwer, 1999)
- Berry & Keating, "The Riemann Zeros and Eigenvalue Asymptotics," *SIAM Review* 41(2), 1999

### The Operator

**Classical Hamiltonian:** H = xp (position times momentum on the half-line)

**Symmetrized quantum operator (Weyl quantization):**
```
H_BK = (1/2)(x̂p̂ + p̂x̂) = -i(x d/dx + 1/2)
```
This is the generator of dilations on R_+, acting on L²(R_+, dx).

### The Idea

The classical trajectories of H = xp are hyperbolas xp = E in phase space. The key observation is that the Bohr-Sommerfeld quantization of this system, after regularization, reproduces the smooth (Weyl) part of the zero counting function:

```
N̄(E) = (E/2π) log(E/2πe) + 7/8
```

This matches the leading asymptotics of the Riemann-von Mangoldt formula for the number of zeros up to height E.

### Regularization

The classical system has unbounded orbits — the phase space is non-compact, so there are no discrete eigenvalues. Berry-Keating proposed truncating phase space:
- |x| ≥ ℓ_x (position cutoff)
- |p| ≥ ℓ_p (momentum cutoff)
with ℓ_x · ℓ_p = 2πℏ (minimum uncertainty cell).

A particle starts at (ℓ_x, E/ℓ_x), follows the hyperbola, and ends at (E/ℓ_p, ℓ_p) after time T = log(E/(ℓ_x ℓ_p)). The semiclassical counting function from this truncation gives exactly N̄(E) including the 7/8 constant.

### What Is Proven
- The smooth counting function N̄(E) is exactly reproduced semiclassically.
- The operator -i(x d/dx + 1/2) on L²(R_+, dx) has purely continuous spectrum (proven by Endres & Steiner, 2009), so it cannot directly yield discrete zeros.
- On compact quantum graphs, all self-adjoint extensions of H_BK have discrete spectra, but the Weyl asymptotics are wrong — they cannot match zeta zeros (Endres & Steiner, 2009).

### Where/Why It Fails
1. **No discrete spectrum on the natural space.** The operator on L²(R_+) has continuous spectrum only.
2. **Regularization is ad hoc.** Truncating phase space reproduces the average density but not individual zeros.
3. **The oscillatory part d_osc(E)** (encoding individual zeros via prime numbers) is not captured.
4. **No natural boundary condition** has been found that makes the operator both self-adjoint and yields the correct discrete spectrum.
5. **The integer dilation condition** f(nx) = f(x) is suggestive but has not led to a rigorous construction.

### Significance
Despite failing as a concrete operator, the Berry-Keating conjecture provides the strongest physical intuition: any successful Hilbert-Polya operator should have H = xp as its classical limit. This constraint has guided all subsequent work.

---

## 2. Berry-Keating Generalizations: H = x(p + ℓ²_p/p) and Relatives

### Key References
- Sierra & Rodriguez-Laguna, "H = xp Model Revisited and the Riemann Zeros," *Phys. Rev. Lett.* 106, 200201 (2011)
- Berry & Keating, "A compact Hamiltonian with the same asymptotic mean spectral density as the Riemann zeros," *J. Phys. A* 44, 285203 (2011)

### The Operators

Several groups independently found that adding a confining term to xp produces closed orbits and discrete spectra:

**Sierra-Rodriguez-Laguna (2011):**
```
H = x(p + ℓ²_p / p)
```
This has closed periodic orbits. The quantized version's smooth spectral density matches N̄(E) for the Riemann zeros.

**Berry-Keating (2011):**
```
H = (x + ℓ_x)(p + ℓ_p)    [compact version]
```
Expanding: H = xp + ℓ_p x + ℓ_x p + ℓ_x ℓ_p. The phase-space area enclosed by the energy-E orbit gives the correct counting function when ℓ_x ℓ_p = 2πℏ.

**The general family:**
```
H = U(x)p + V(x)/p
```
was studied for various U, V. For U(x) = x, V(x) = ℓ²_p x, one recovers the Sierra model. The dynamics describe a massive particle in a relativistic spacetime whose metric is constructed from U and V.

### What Is Proven
- The smooth counting function N̄(E) is reproduced for several choices.
- The models have compact classical orbits and therefore discrete quantum spectra.
- For H = (x + ℓ_x)(p + ℓ_p), the Weyl asymptotics match exactly.

### Where It Fails
- Only the **average** zero density is reproduced, not individual zeros.
- The **oscillatory corrections** (which encode prime numbers) are absent.
- No choice of U(x), V(x) has been found that reproduces the exact spectrum.
- The models are classically integrable, creating tension with GUE statistics (which indicate chaotic dynamics).

---

## 3. Sierra-Townsend Landau Level Model (2007-2008)

### Key References
- Sierra, "H = xp with interaction and the Riemann zeros," *Nucl. Phys. B* 776, 327 (2007), arXiv:math-ph/0702034
- Sierra & Townsend, "Landau levels and Riemann zeros," *Phys. Rev. Lett.* 101, 110201 (2008)

### The Operator

**Sierra (2007) — xp with interaction:**

Starting from H = xp, Sierra adds a non-local interaction depending on two wave functions (potentials). The model is formulated via a Jost-like function:

```
H = xp + V_interaction(x)
```

The interaction is chosen so that the Jost function φ(k) is analytic in the upper half-plane, with:
- Real zeros on the real axis → bound states
- Complex zeros below the real axis → resonances

For specific potentials, the resonances converge asymptotically toward the average positions of the Riemann zeros.

**Key result:** A linear superposition of potentials, obtained by the action of integer dilations, yields a Jost function whose real part vanishes at the Riemann zeros and whose imaginary part approximates Im ζ(1/2 + it).

**Sierra-Townsend (2008) — Landau levels:**
```
H = (1/2μ)[p_x² + (p_y + eBx/c)²] + eλxy
```
A charged particle on a 2D surface in:
- A uniform perpendicular magnetic field B
- An electric potential V = λxy

In the lowest Landau level limit, this reduces to the Connes "absorption spectrum" model. The smooth counting function N̄(E) is recovered from the Landau level structure.

### Space
- 2D plane (x, y) for the Landau level model
- Half-line R_+ for the 1D Jost function model

### What Is Proven
- The average zero density is recovered in the lowest Landau level.
- The Jost function construction produces bound states at positions matching average Riemann zeros.
- The x ↔ p exchange symmetry of H = xp is exploited to make Jost solutions proportional to ζ(s).

### Where It Fails
- Sierra explicitly states: "We cannot exclude the existence of zeros outside the critical line."
- Higher Landau level contributions to the oscillatory terms are only speculatively connected to fluctuation corrections.
- The Jost function approach requires inputting the Riemann-Siegel formula.
- The model suggests zeros as resonances embedded in a continuum, but proving discreteness requires additional assumptions.

---

## 4. Bender-Brody-Muller PT-Symmetric Hamiltonian (2017)

### Key References
- Bender, Brody & Muller, "Hamiltonian for the Zeros of the Riemann Zeta Function," *Phys. Rev. Lett.* 118, 130201 (2017), arXiv:1608.03679
- Bellissard, "Comment on 'Hamiltonian for the Zeros of the Riemann Zeta Function,'" arXiv:1704.02644 (2017)
- Bender, Brody & Muller, reply, arXiv:1705.06767 (2017)

### The Operator

```
Ĥ = (1 - e^{-iπ̂})^{-1} · (x̂p̂ + p̂x̂) · (1 - e^{-ip̂})
```

More precisely, it is a similarity transformation of the Berry-Keating operator:
```
Ĥ = Ŝ⁻¹ · H_BK · Ŝ
```
where Ŝ = (1 - e^{-ip̂}).

### Properties
- **Classical limit:** H → 2xp (consistent with Berry-Keating)
- **Not Hermitian** in the conventional sense
- **iĤ is PT-symmetric** with broken PT symmetry
- PT symmetry allows the possibility that all eigenvalues are real

### The Claim
If the eigenfunctions obey the boundary condition ψ(0) = 0, then the eigenvalues correspond to the nontrivial zeros of the Riemann zeta function. The authors present a heuristic argument for constructing a metric operator (inner product) under which Ĥ becomes Hermitian. If this can be made rigorous, RH follows.

### Bellissard's Criticism (2017)
Bellissard showed that the construction has a fatal domain problem:
1. The similarity transformation Ŝ = (1 - e^{-ip̂}) involves the **momentum operator** p̂.
2. On L²(R_+) with the boundary condition ψ(0) = 0, the momentum operator p̂ = -i d/dx is **not self-adjoint** (it has deficiency indices (1,0) or (0,1) on the half-line).
3. Therefore e^{-ip̂} is not a well-defined unitary operator on this domain.
4. The similarity transformation Ŝ is **not well-defined** on the domain where H_BK is self-adjoint and the boundary condition applies.
5. The resulting Hamiltonian Ĥ is therefore not a legitimate operator on the intended Hilbert space.

### Authors' Response
Bender, Brody & Muller responded that Bellissard's objections had already been addressed in their paper and do not affect the conclusions. However, the mathematical community has not found this response convincing.

### Current Status
The BBM approach remains **unresolved but widely regarded as incomplete**. The fundamental domain issue identified by Bellissard has not been rigorously overcome. The paper is interesting as a suggestion but falls far short of a proof. Subsequent work (e.g., the 2023 paper arXiv:2309.00405) has proposed alternative similarity transformations using the **number operator** and **position operator** (both self-adjoint on L²[0,∞)) to avoid the momentum domain problem.

---

## 5. Connes' Noncommutative Geometry Approach (1996-2024)

### Key References
- Connes, "Trace formula in noncommutative geometry and the zeros of the Riemann zeta function," *Selecta Math.* 5, 29-106 (1999), arXiv:math/9811068
- Connes, "An essay on the Riemann Hypothesis," arXiv:1509.05576 (2015/2019)
- Connes, Consani & Moscovici, "Zeta zeros and prolate wave operators," *Ann. Funct. Anal.* (2024), arXiv:2310.18423
- Connes & Consani, "Zeta spectral triples," arXiv:2511.22755 (2025)

### The Framework

Connes does not propose a single Hamiltonian. Instead, he constructs a geometric framework using noncommutative geometry:

**The space:** The adele class space X = A_Q / Q*, where A_Q is the ring of adeles of Q and Q* acts by multiplication. This is a noncommutative space (the quotient is "bad" in classical geometry).

**The action:** The multiplicative group R*_+ acts on X by scaling. This action generates a one-parameter flow, and the "frequencies" of this flow are the zeros of zeta.

**Spectral interpretation:** The nontrivial zeros appear as an **absorption spectrum** — they are missing frequencies in the continuous spectrum, rather than discrete eigenvalues. This is dual to the Hilbert-Polya picture: instead of eigenvalues of an operator, the zeros are points where a continuous family of representations fails to appear.

### The Trace Formula

Connes proves that the Weil explicit formula (relating zeros of zeta to primes) can be rewritten as a **trace formula** on the adele class space, analogous to the Selberg trace formula for hyperbolic surfaces:

```
∑_ρ ĥ(ρ) = ĥ(0) + ĥ(1) - ∑_p ∑_{m=1}^∞ (log p / p^{m/2}) h(m log p) + (distributional terms)
```

The left side sums over zeta zeros ρ; the right side sums over primes p.

### What Is Proven
- The trace formula is **equivalent to the Riemann Hypothesis** for all L-functions with Grossencharakter: RH holds if and only if the trace formula is valid in the global case.
- The absorption spectrum interpretation resolves a sign discrepancy in the Gutzwiller trace formula approach.
- The framework unifies the Riemann-Weil explicit formula with the Selberg trace formula in a single geometric picture.

### Where It Stalls
1. **The positivity problem.** Proving the trace formula amounts to proving a positivity condition (Weil positivity) for a certain distribution. This positivity is equivalent to RH but has not been established.
2. **No explicit operator with computable eigenvalues.** The framework is geometric/algebraic, not providing a concrete Hamiltonian one could diagonalize.
3. **The "missing piece" is essentially as hard as RH itself.** Connes has reformulated RH as a statement in noncommutative geometry but has not reduced its difficulty.

### Recent Progress: Prolate Wave Operators (2024)

In the latest work, Connes (with Consani and Moscovici) introduces a **semilocal prolate wave operator**:
- In the archimedean case: the prolate operator = (scaling operator)² + (grading of orthogonal polynomials)
- The **positive spectrum** of this operator realizes the low-lying zeta zeros
- The **negative spectrum** corresponds to the Sonin space (governing ultraviolet behavior)
- The operator is related to the **metaplectic representation** of the double cover of SL(2,R)

This represents the most concrete operator construction to come out of Connes' program, though it still does not constitute a proof.

### Zeta Spectral Triples (2025)

The most recent paper constructs **self-adjoint operators D(λ,N)** as rank-one perturbations of a spectral triple associated with scaling operators on [λ⁻¹, λ], involving Euler products over primes p ≤ N. These produce self-adjoint operators whose spectra coincide numerically with the lowest nontrivial zeros of ζ(1/2 + is). This is the closest any construction has come to an explicit Hilbert-Polya operator, but rigorous proofs of spectral coincidence in the limit are still missing.

---

## 6. Wu-Sprung Potential (1993)

### Key References
- Wu & Sprung, "Riemann zeta zeros and the inverse problem," *Phys. Rev. E* 48, 2595 (1993)

### The Operator

A standard 1D Schrodinger operator:
```
H = -d²/dx² + V_WS(x)
```
on L²([0,∞)) with Dirichlet boundary condition ψ(0) = 0.

The potential V_WS(x) is determined **inversely** from the requirement that the eigenvalues match the Riemann zeros. Using the WKB (semiclassical) inversion formula:

```
x(V) = (1/π)[√(V - V₀) ln(V₀/2π) + √V · ln((√V + √(V-V₀))/(√V - √(V-V₀)))]
```

For large x:
```
x(V) ≈ (√V / π) ln(2V / πe²)
```

### What Is Proven
- The smooth semiclassical density-of-states matches N̄(E) for the first ~500 zeros.
- The potential has **fractal structure** with dimension d ≈ 1.5.
- Multiple independent reconstruction methods (WKB inversion, Gel'fand-Levitan, Marchenko) yield the same potential, confirming uniqueness in 1D.
- The fractal/multi-fractal character has been confirmed by independent groups (Schumayer, van Zyl, Hutchinson, etc.).

### Where It Fails
- **Fundamental contradiction with GUE:** The Wu-Sprung model is a 1D integrable system. According to Berry-Tabor/BGS conjectures, integrable systems should show Poisson statistics, not GUE. Yet the potential is specifically designed to produce GUE-distributed eigenvalues. This works because the potential is fractal (non-smooth), effectively bypassing the integrability argument — but this loophole is poorly understood.
- **The potential is not smooth:** It is a fractal function, with no closed-form expression. This makes it physically unrealizable and mathematically unwieldy.
- **Circular reasoning risk:** The potential is constructed from the zeros, so proving RH from this operator would require proving something about the fractal structure of V_WS that is essentially equivalent to RH.
- **Only the smooth part determines the potential.** The oscillatory corrections require least-squares fitting to individual zeros.

### Significance
Wu-Sprung is important as an **existence proof**: a self-adjoint operator with the Riemann zeros as eigenvalues exists (by inverse spectral theory). The question is whether it has additional structure (symmetry, number-theoretic content) that would constitute a proof of RH.

---

## 7. Bhaduri-Ghosh-Vyas Inverted Oscillator (1995-1997)

### Key References
- Bhaduri, Khare & Law, *Phys. Rev. Lett.* 74, 4963 (1995)

### The Operator

**1D version:**
```
H = -d²/dx² - (1/2)mω²x²    (inverted harmonic oscillator)
```
on L²([0,∞)) with boundary condition Φ(0) = 0.

**2D version:** Extension with a perpendicular parabolic confining potential and inverted oscillator in the y-direction, incorporating Landau-level-like structure.

### What Is Proven
- The oscillating scattering phase shift δ(t) of the inverted oscillator equals exactly the phase of ζ(1/2 + it).
- The connection is via repulsive Coulomb scattering mapping.
- The Argand diagram of the scattering amplitude resembles that of the zeta function.

### Where It Fails
- Only reproduces the oscillatory part of the zeta phase, not the full zero structure.
- The 2D model's fluctuation terms are only supported at order-of-magnitude level.
- The boundary condition at the origin is artificial with no number-theoretic motivation.

---

## 8. Chadan-Khuri Model (1993-2002)

### Key References
- Khuri, "A model for the Riemann zeros," *Math. Phys. Anal. Geom.* 5, 1 (2002)

### The Operator

```
H = -d²/dr² - ℓ(ℓ+1)/r² + (1/r²) f_CM(r)
```
A radially symmetric 3D Hamiltonian, restricted to r ∈ [0, e^{-4π/3}].

### What Is Proven
- The coupling constant spectrum approximately coincides with nontrivial zeros.
- A 3D potential is derived whose s-wave scattering amplitude has zeros on the critical line.

### Where It Fails
- Only "approximate" agreement with zeros.
- The finite interval restriction is artificial.
- Logarithmic singularities at r = 0 lack physical interpretation.

---

## 9. Crehan's Existence Theorem (1995)

### Key Reference
- Crehan, "Chaotic spectra of classically integrable systems," *J. Phys. A* 28, 6389 (1995)

### Result
For **any** bounded sequence of real numbers {E_n}, there exist infinitely many classically integrable nonlinear oscillator Hamiltonians whose quantum spectra match {E_n}. Applied to the Riemann zeros as a special case.

### Significance
This is a pure existence result — non-constructive, giving no explicit operator. It shows that the existence of a Hilbert-Polya operator is trivially guaranteed (in a non-useful sense). The real question is finding an operator with **additional structure** (number-theoretic content, spectral symmetries).

---

## 10. Supersymmetric Quantum Mechanics Approaches (2019-2025)

### Key References
- Ramos, Arias de Saavedra & Falceto, "Supersymmetry and the Riemann zeros on the critical line," *Phys. Lett. B* 795, 418 (2019)
- Kar, "Supersymmetric quantum mechanical system for locating the Riemann zeros," *Eur. Phys. J. Plus* 138, 487 (2023)
- De Angelis, "A Supersymmetric Quantum mechanical model and the spectral embedding conjecture for the Riemann zeros," (2025)

### The Idea

Rather than constructing an operator whose **entire** spectrum is the set of Riemann zeros, the SUSY approach embeds the zeros within a larger spectrum:

A family of **supersymmetric quantum mechanical Hamiltonians** on the half-line [0,∞) features:
- A confining **logarithmic potential** (ensuring discrete spectrum)
- A **scale-invariant conformal core** (reflecting dilation symmetry of xp)
- **Symmetry-breaking perturbations** (encoding arithmetic information)

The Witten superpotential W(x) is chosen so that the SUSY ground state energy E₀ = 0 occurs precisely when the spectral parameter s is a nontrivial zero of ζ(s).

### What Is Proven
- SUSY partner Hamiltonians H± can be constructed whose vanishing ground-state energy locates zeros on the critical line.
- The spectral embedding conjecture (zeros embedded in a continuous family of spectra) provides a new framework.

### Where It Fails/Stalls
- The approach reframes the Hilbert-Polya problem as a **spectral selection principle** — finding which parameter values yield E₀ = 0 is equivalent to finding zeros of ζ.
- No proof that the spectral parameter space is restricted to Re(s) = 1/2.
- The models have not been shown to exclude off-line zeros.

---

## 11. The 2023 Riemann Operator Construction (arXiv:2408.15135)

### Key Reference
- "Nontrivial Riemann Zeros as Spectrum," arXiv:2408.15135 (2024)

### The Operator

```
R̂ = -D̂ - iμ(T̂)
```
where:
- D̂ = Berry-Keating Hamiltonian = -i(x d/dx + 1/2)
- T̂ = Bessel operator (self-adjoint)
- μ(T̂) = T̂ tanh(T̂/2) - I

**Domain:** D(R̂) = {ψ ∈ D(D̂) ∩ D(T̂) | ψ(0) = 0} in L²([0,∞), dx)

### Eigenstates
```
|Ψ_λ⟩ = ∫₀^∞ t^{λ-1} ω(t) |t⟩ dt
```
where ω(t) = t e^t / (1 + e^t)² is related to the Fermi-Dirac distribution.

### Self-Adjoint Version
R̂ itself is **non-self-adjoint**. A self-adjoint Hilbert-Polya operator is constructed via similarity:
```
ĥ = Ŵ^{1/2} · R̂_{S_ζ} · Ŵ^{-1/2}
```
where Ŵ is a positive semidefinite intertwining operator satisfying V̂R̂ = R̂†V̂.

Essential self-adjointness is established via deficiency index analysis: dim ker(ĥ† ∓ i) = 0.

### Key Theorem (5.1)
If the intertwining operator satisfies Ŵ ≥ 0, then all simple zeros satisfy Re(ρ) = 1/2. This reduces RH to proving positivity of a specific operator.

### Non-Uniqueness
The construction is explicitly non-unique (Remark 3.3): the Bessel operator may be replaced by any self-adjoint operator whose generalized eigenfunctions satisfy completeness and homogeneity. The weight function ω(t) can also be replaced, as long as the Mellin transform includes zeta zeros and eigenfunctions remain L².

---

## The Inverse Spectral Problem

### What Constraints Does Self-Adjointness + Matching Zeros Impose?

Given the Riemann zeros {γ_n} as a prescribed spectrum, the inverse spectral problem asks: what self-adjoint operators have this spectrum?

**Key constraints that can be attacked independently (from arXiv:2408.15135):**
- **(A1) Self-adjointness:** The operator must be self-adjoint (or admit only real eigenvalues via PT symmetry or similar).
- **(A2) Spectral matching:** Eigenvalues must coincide with {γ_n} (explicit formula matching).
- **(A3) Functional equation symmetry:** The operator must implement an involution reflecting the functional equation ξ(s) = ξ(1-s).
- **(A4) Determinant identity:** ξ(s) must arise as the spectral determinant (regularized Fredholm determinant) of the operator.

### Non-Uniqueness

The inverse spectral problem is **highly non-unique** in general:
- In 1D, specifying a single spectrum determines the potential only up to the choice of norming constants (Borg-Marchenko theory). Two spectra (e.g., Dirichlet and Neumann) uniquely determine the potential.
- Crehan's theorem guarantees infinitely many classically integrable Hamiltonians for any spectrum.
- The Wu-Sprung construction shows multiple reconstruction methods yield the same 1D potential (uniqueness in the 1D Schrodinger class), but this is one specific class among many.

### What Narrows the Space?

**GUE statistics (Montgomery-Odlyzko):** The pair correlation of zeta zeros matches GUE. This implies:
- The operator should **break time-reversal symmetry** (GUE rather than GOE or GSE).
- The classical limit should be **chaotic** (not integrable — though Wu-Sprung shows fractal potentials can evade this).
- Short-range correlations follow universal GUE predictions.
- Long-range correlations deviate from GUE universality beyond scale ~ln(E/2π)/ln 2, reflecting system-specific (prime-number) structure.

**Functional equation symmetry:** The involution s ↔ 1-s must be implemented as an operator symmetry. In Berry-Keating language, this corresponds to x ↔ p exchange symmetry.

**Spectral determinant = ξ:** The completed zeta function ξ(s) should be the Fredholm/spectral determinant of the operator, encoding both the spectrum and its multiplicities.

**Prime structure in the trace formula:** The periodic orbit sum (or its analogue) must reproduce the prime numbers, via a trace formula analogous to Selberg's.

---

## The Parameterized Family H = f(x)p + g(x)

### What's Known

The general family:
```
H = U(x)p + V(x)/p
```
was studied systematically. This describes a massive relativistic particle in a spacetime whose metric is determined by U and V.

**Known cases:**
| U(x) | V(x) | Name | Result |
|-------|-------|------|--------|
| x | 0 | Berry-Keating | Correct N̄(E), continuous spectrum |
| x | ℓ²_p x | Sierra 2011 | Closed orbits, correct N̄(E) |
| x + ℓ_x | 0 (with p + ℓ_p) | Berry-Keating compact | Correct N̄(E) + 7/8 |
| x + 1/x | (with p + 1/p) | Sierra-Rodriguez-Laguna | Correct smooth asymptotics |

### Has Anyone Searched Computationally?

There is no published systematic computational search over the (U, V) parameter space for operators matching individual Riemann zeros. The literature has focused on analytic/semiclassical matching of the smooth counting function. A computational search for exact eigenvalue matching would require:
1. Choosing a parametric form for U, V
2. Solving the quantum eigenvalue problem numerically
3. Comparing eigenvalues to known zeros
4. Optimizing over parameter space

This appears to be an **open direction** that has not been pursued in the literature, likely because:
- The smooth part is easy to match (many models do it)
- The oscillatory part (individual zeros) requires encoding prime numbers, which no simple parametric family captures
- Any match found computationally would need a proof to be meaningful

---

## Summary Table

| # | Candidate | Year | Operator | Status |
|---|-----------|------|----------|--------|
| 1 | Berry-Keating | 1999 | H_BK = -i(x d/dx + 1/2) | Continuous spectrum on L²(R_+); correct smooth asymptotics only |
| 2 | Berry-Keating generalizations | 2011 | H = x(p + ℓ²/p), etc. | Discrete spectra, correct average density, no individual zeros |
| 3 | Sierra-Townsend | 2007-08 | xp + interaction; Landau levels | Zeros as resonances; cannot exclude off-line zeros |
| 4 | Bender-Brody-Muller | 2017 | Ŝ⁻¹ H_BK Ŝ, PT-symmetric | Domain problems (Bellissard); unresolved |
| 5 | Connes NCG | 1996-2024 | Adele class space trace formula; prolate wave operators | RH ⟺ trace formula validity ⟺ Weil positivity; not proven |
| 6 | Wu-Sprung | 1993 | -d²/dx² + V_fractal(x) | Existence proof; fractal potential; circular |
| 7 | Bhaduri et al. | 1995 | Inverted harmonic oscillator | Phase matches; not full spectrum |
| 8 | Chadan-Khuri | 1993-2002 | Radial Schrodinger, coupling constant | Approximate only |
| 9 | Crehan | 1995 | (existence theorem) | Non-constructive |
| 10 | SUSY QM | 2019-2025 | Partner Hamiltonians, log potential | Spectral embedding; no exclusion of off-line zeros |
| 11 | Riemann operator (2024) | 2024 | R̂ = -D̂ - iμ(T̂) | RH ⟺ positivity of Ŵ; most explicit construction to date |
| 12 | Connes zeta spectral triples | 2025 | D(λ,N) rank-one perturbations | Numerical agreement; rigorous limit unproven |

---

## Key Open Questions

1. **Can the oscillatory part of the counting function be recovered from any concrete operator?** All current models match only the smooth Weyl term. The oscillatory corrections encode the primes via a trace formula, and no proposed Hamiltonian naturally produces them.

2. **Is the Hilbert-Polya operator unique (up to unitary equivalence) once all constraints are imposed?** The inverse spectral problem is non-unique in general, but the combined constraints of self-adjointness + GUE statistics + functional equation symmetry + spectral determinant = ξ may determine it uniquely.

3. **What is the correct Hilbert space?** Candidates include L²(R_+), L²(adele classes), spaces of automorphic forms, and Hilbert spaces of entire functions. The choice of space is inextricable from the choice of operator.

4. **Can the positivity conditions (Connes' Weil positivity, or Ŵ ≥ 0 in the 2024 construction) be proven?** Both known reformulations reduce RH to positivity of a specific mathematical object.

5. **Is there a natural physical system whose quantum Hamiltonian is the Hilbert-Polya operator?** The GUE statistics suggest a quantum chaotic system without time-reversal symmetry. Berry has speculated it might involve a charged particle in a magnetic field — consistent with Sierra-Townsend's Landau level approach.

---

*Last updated: 2026-03-30*
*Part of Month 3-4 research phase*
