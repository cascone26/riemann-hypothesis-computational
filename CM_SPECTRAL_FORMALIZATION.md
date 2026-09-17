# Machine-Verified Formalization of Connes–Moscovici's Spectral Realization
*(Weil-form → ζ-spectrum construction), as of 2026-08-30 (cy376–390)*

**Status: axiom-clean and real — a Lean 4 + Mathlib formalization of the *conditional* spectral-realization
lemma from Connes–Moscovici (arXiv:2511.22755, "ccm"). It is HANDLE-tier (compiled, audited). It is NOT a
proof of RH: the construction is conditional on CM's own open "even-simple" hypothesis (their §1385 missing
step #1), which is superexponentially hard and numerically inaccessible (see §D).**

> **What this is:** the first machine-checked formalization of CM's construction that turns the minimal
> eigenvector `ξ` of the truncated Weil quadratic form into a self-adjoint operator `D''` whose spectrum is
> the zero set of `ξ̂` — assembled abstractly, then instantiated on the *actual* Weil matrix with every
> structural hypothesis discharged. **What it is not:** a closure of CM's program. The lone open hypothesis
> (even-simple) is exactly the analytic wall the field itself has not crossed.

> **This is the finite core of Connes–van Suijlekom Theorem 6.1 (verified cy430).** Connes 2026
> (arXiv:2602.04022, "The Riemann Hypothesis: Past, Present and a Letter Through Time", §6) states **Theorem
> 6.1** (joint with W. van Suijlekom, ref [32] "Quadratic Forms, Real Zeros and Echoes of the Spectral
> Action"): *for an even distribution whose Weil-type quadratic form is lower-bounded self-adjoint with a*
> ***simple isolated minimal eigenvalue and even eigenfunction η, all zeros of η̂ lie on the real line.*** Its
> proof is, in Connes's own words, *"the special form of the matrix of the quadratic form in the trigonometric
> orthonormal basis, the construction for finite matrices of that special form of a selfadjoint operator, and
> Hurwitz to pass to the limit ... approximate the zeros of η̂ by the spectrum of a rank one perturbation of
> the periodic Dirac operator."* That finite construction is **exactly** our master theorem (τ-form → self-
> adjoint `D'' = D − |Dξ⟩⟨η|`, spectrum = zeros of `ξ̂`) plus our compression chain (simple isolated minimum
> ⟹ PSD form with 1-dim radical `⟨ξ⟩` ⟹ PD on quotient ⟹ real spectrum). So the forge machine-verifies the
> **finite core of a named 2026 Connes–van Suijlekom theorem** (modulo the basis-existence assembly step). And
> Connes 2026 §6.6 ("Remaining steps") states the two open inputs verbatim — *"smallest eigenvalue of QW_λ is
> simple with even eigenvector"* and *"k_λ is a sufficiently good approximation of θ_x"* — independently
> confirming this document's honest scope; §6.4 confirms *"RH implies QW_λ is strictly positive and its radical
> is {0}"*, which is exactly the compression's radical-trivial ⟺ PD ⟺ real-spectrum mechanism.

Location: Olympus `~/rhforge/forge/` (namespace `CMKeyLemma`). Axiom gate: `{propext, Classical.choice,
Quot.sound}` only — re-audited GREEN 2026-08-30 (`AXIOM AUDIT: checked=1 clean=1 bad=0` on
`WeilSpectralRealization`, `DetDprimeDiagonal`, `QDfromForm`, and every module below at bank time).

> **Exact lemma-level correspondence (verified cy431) to Connes–van Suijlekom §5** (arXiv:2511.23257,
> "Quadratic Forms, Real Zeros and Echoes of the Spectral Action", 28 Nov 2025 — sibling of the ccm paper):
> our τ-form = their `q_ij=(b_i−b_j)/(λ_i−λ_j)`, which their §4.3 identifies as **literally the second Gâteaux
> derivative (Hessian) of the spectral action**. `QDfromForm`+`GammaCommute` = their **Lemma 5.1**;
> `CMKeyLemmaI`+`CMKeyLemmaII` = their **Lemma 5.2**; `DetDprimeDiagonal` = their **Lemma 5.3**; our compression
> chain (`posDef_compress`+capstone+`real_spectrum_blocktri`) = their **Theorem 5.6(i) proof** step-for-step
> (PSD 1-dim radical `⟨ξ⟩` → self-adjoint `D''` on the separated quotient → real spectral theorem → triangular
> basis `{ξ, w_j}` → `Det(D'−s)=−s∏(λ_j−s)`, all real); the master theorem's L-agnostic generality = their
> **§5.1 Prop 5.10** (general simple-spectrum `D`). The deferred basis-existence step is exactly their
> triangular basis `{ξ, w_j}` (`w_j` = lifts of the `D''`-eigenvectors) — construction blueprint now explicit.

---

## A. The master theorem (`WeilSpectralRealization.weil_spectral_realization`)

Let `D = diagonal c` be the scaling operator (`c` = the distinct eigenvalues, CM's `c_j = j`) and `T` the
Weil quadratic form `QW_λ^N` in CM's τ-form. Given the grading symmetry (`a` even, `b`,`c` odd under the
reversal involution `r`) and **even-simple** as the sole external input (`T·ξ = 0`, `ξ` even, `⟨η,ξ⟩ = 1`,
`η ≡ 1`), for `D' = D − |Dξ⟩⟨η|` and any `s` off the spectrum (`c_j − s ≠ 0`):

> **(1)** `(T·D')ᵀ = T·D'`  — `D''` is self-adjoint in the `T`-inner product.
> **(2)** `det(D'−s·1) = −s · det(D−s·1) · Σ_j (c_j−s)⁻¹ ξ_j` — with `c_j=j`, `η≡1` the sum is `ξ̂(s)`, so
> the eigenvalues of `D''` are exactly the **zeros of `ξ̂`**.

All of {real-symmetry, the QD commutator, `⟨β,ξ⟩=0`, `D` symmetric} are **discharged internally** from the
τ-form and grading. Only even-simple is assumed.

## B. The brick inventory (14 modules, all axiom-clean)

**Abstract key lemma (CM Lemma "key" / §850–938):**
| Module | Statement |
|---|---|
| `CMKeyLemmaI` | `T·Dξ = −β` from the commutator + `Tξ=0, ⟨η,ξ⟩=1, ⟨β,ξ⟩=0` (part i) |
| `CMKeyLemmaII` | `(T·D')ᵀ = T·D'` — self-adjointness of `D''` (part ii); consumes I |
| `DetAddVecMulVec` | matrix determinant lemma `det(A+|u⟩⟨v|) = detA·(1+⟨v,A⁻¹u⟩)` (engine) |
| `ResolventDXi` | `(D−s)⁻¹·Dξ = ξ + s(D−s)⁻¹ξ` (resolvent identity) |
| `DetDprimeCore` | `det(D'−s) = −s·det(D−s)·⟨η,(D−s)⁻¹ξ⟩` (part iii core); composes the two engines |
| `DiagResolventForm` | `⟨η,(diag g)⁻¹ξ⟩ = Σ_j η_j(g_j)⁻¹ξ_j` (ξ̂ connector) |
| `DetDprimeDiagonal` | `det(D'−s) = −s·det(D−s)·Σ_j η_j(c_j−s)⁻¹ξ_j` (part iii capstone) |
| `DprimeKernel` | `D'ξ = 0` (kernel fact founding the `D''` quotient / det0) |

**Concrete instantiation (CM Lemma "basics" / eq. `form`, §819):**
| Module | Statement |
|---|---|
| `QDfromForm` | the τ-form ⟹ the QD commutator `DT−TD = \|β⟩⟨η\| − \|η⟩⟨β\|` (concrete→abstract bridge) |
| `GammaCommute` | the τ-form (a even, b,c odd) ⟹ `T` commutes with γ (`T.submatrix r r = T`) |
| `OddEvenOrtho` | `b` odd, `ξ` even ⟹ `⟨b,ξ⟩ = 0` (discharges `hβ`) |
| `TSymmOfForm` | the τ-form ⟹ `Tᵀ = T` (real-symmetry) |

**Master + a sibling:**
| Module | Statement |
|---|---|
| `WeilSpectralRealization` | the master theorem above (composes all of the above) |
| `SoninArchWeil` | (cy370, independent) archimedean Weil positivity from the Sonin trace bound: `Tr(V𝔖Vᴴ)≤W_∞ ⟹ W_∞≥0` |

## C. Composition graph

```
  QDfromForm ─┐                              (concrete τ-form ⟹ QD commutator)
  TSymmOfForm ─┤                              (τ-form ⟹ Tᵀ=T)
  OddEvenOrtho ┤──► CMKeyLemmaII ──┐          (self-adjointness, part ii)
  CMKeyLemmaI ─┘        ▲          │
                        (part i)   ├──► weil_spectral_realization  [MASTER]
  DetAddVecMulVec ─┐               │
  ResolventDXi ────┴─► DetDprimeCore ─► DetDprimeDiagonal ─┘   (det realization, part iii)
  DiagResolventForm ──────────────────────┘
  DprimeKernel  (D'ξ=0, foundation of det0 — the D'→D'' passage, not yet formalized)
```

## D. The one open hypothesis (honest status)

The master assumes **even-simple**: the smallest eigenvalue of `QW_λ` is simple with an *even* eigenvector.
This is **CM's own stated open "missing step" #1** (ccm §1385) — proven by CM for the auxiliary *prolate*
operator `PW_λ`, but *not* for the Weil form `QW_λ`. Our independent work confirms it is **numerically
inaccessible**:
- the true `ε_λ` is superexponentially small (CM report `log₁₀ ε_λ ≈ −530`);
- a faithful high-precision build of `QW_λ` (`operators/qw_lambda_faithful.py`, sign-verified against CM eq.
  `bombtest`, grading-commuting to `‖Q_eo‖≈5e-16`) floors out — even the even-vs-odd sector-minimum
  comparison is float64 noise whose verdict flips with `λ` (cy383).

So the transfer `PW_λ` simple-even ⟹ `QW_λ` simple-even needs an **analytic** mechanism (CM's Co-zeta
trace-formula bridge), which is external-math-scale, not autonomous-closable. This formalization isolates
that hypothesis as the *lone* input to an otherwise-complete, machine-verified construction.

## E. Not yet formalized (deferred)
- `det0`: `Det(D'−s) = −s·Det(D''−s)` — the explicit `D'→D''` quotient-operator determinant (needs Mathlib
  quotient-determinant infrastructure; heavy). `DprimeKernel` (`D'ξ=0`) is its foundation.
- the §835 normalization derivation (`even-simple ⟹ ⟨η,ξ⟩ ≠ 0`), which would derive the normalization from
  the 1-dim kernel rather than assume it (marginal — `⟨η,ξ⟩=1` is a standard eigenvector scaling).

_Source: THE_BRIDGE.md cy376–390; forge `~/rhforge/forge/` on Olympus; CM = arXiv:2511.22755._

---

## Scope: the master theorem proves GRH, not just RH (verified cy415)

`weil_spectral_realization` was checked to have **no ζ-specific hypothesis** — its signature takes only abstract
τ-form data `(c a b ξ : n → ℝ)`, a reversal involution `r`, and the matrix `T` with the form/grading/even-simple
hypotheses. It is pure linear algebra on *any* symmetric τ-form matrix. The Weil explicit-formula structure is
universal across the Selberg/Dirichlet class — only the entries `a,b,c` (archimedean shift + prime weights)
change per L-function. Therefore the SAME theorem gives, for any Dirichlet L(s,χ):

> **GRH(χ) ⟸ even-simple(χ)** — machine-verified, axiom-clean, identical proof.

The truncated-Weil numerical pipeline extends to `L(s,χ₃)` with the same structure (truncweil §Future, cy409),
confirming the τ-form is L-function-agnostic. So the deliverable is a **general Selberg/Dirichlet spectral
realization**, conditional on even-simple(χ) — a broader result than "RH alone," with the identical open hinge.

---

## Non-vacuity witness — the master theorem provably applies to a real instance (verified cy417)

`WeilSpectralRealization.weil_spectral_realization` is not vacuously stated: its full hypothesis set
(τ-form + reversal grading + even-simple) is **machine-verified jointly satisfiable** on a concrete,
**non-degenerate** `Fin 3` matrix (`forge/WeilSpectralWitness.lean`, 19th axiom-clean module):

> `c = (−1,0,1)`, `a = (0,2,0)`, `b = (−1,0,1) = c`, `ξ = (1,−1,1)`, `T = [[0,1,1],[1,2,1],[1,1,0]]`,
> `r = Fin.rev`, `s = 1/2`.

Construction: taking `b = c` collapses every off-diagonal τ-entry `(b_i−b_j)/(c_i−c_j)` to `1` (no division
pathology); an even `ξ ∈ ker T` with `⟨1,ξ⟩ = 1` then forces `ξ = (1,−1,1)`, `a = (0,2,0)` (the `α=0` branch —
`α=1` would collapse `T` to the all-ones `J`, whose kernel is `⊥ 1`, making `⟨1,ξ⟩=1` unsatisfiable; the `α=0`
choice is exactly what keeps the normalization hypothesis live). Verified: `Tξ = 0`, `⟨1,ξ⟩ = 1`, `ξ` even,
`a` even, `b` odd, off-diagonal τ-form — all exact. Non-degenerate: `T ≠ 0`, `ξ ≠ 0`, `b ≠ 0`, `c` distinct.

Applying the master theorem to this data typechecks and is axiom-clean:
`#print axioms weil_spectral_realization_nonvacuous = {propext, Classical.choice, Quot.sound}`; the axiom gate
reports GREEN (`checked=20 clean=20 bad=0`, including transitive dependencies). So the deliverable is a theorem
that **provably says something about a real matrix** — non-vacuous, not merely non-contradictory. (This does
*not* touch the open hinge: even-simple for the actual Weil form `QW_λ` remains external/analytic — §D.)

**Module inventory: 19** (was 18 + `WeilSpectralWitness`).

**Conclusion realized end-to-end on the witness (cy418).** Beyond the hypotheses being satisfiable, the master
theorem's *conclusion* was checked to actually fire on this instance. With `D' = diag(c) − |Dξ⟩⟨1|`:
(1) self-adjointness `(T·D')ᵀ = T·D'` holds (proven in Lean; numerics agree); (2) the determinant identity
`det(D'−sI) = −s·det(D−sI)·ξ̂(s)` is exact at every tested `s` (`|diff| ≲ 1e-14`; proven in Lean); and
(3) the **spectral realization** works: `D'ξ = 0` (the `ξ`-kernel), `spec(D') = {i, −i, 0}`, so the quotient
operator `D''` has spectrum `{±i}`, while `ξ̂(s) = Σ_j ξ_j/(c_j−s) = (s²+1)/(s−s³)` has poles at
`c_j = {−1,0,1}` and **zeros exactly at `±i`**. Thus `spec(D'') == zeros(ξ̂)` exactly — the operator's
eigenvalues *are* the zeros of the transform, the CM mechanism demonstrated on a concrete matrix. (The toy
zeros `±i` sit off the real line — a faithful illustration; in the ζ instance `ξ̂`'s zeros are the ζ-zeros
realized the same way.) Facts (1)–(2) are HANDLE (Lean-proven); fact (3)'s spectrum-match is a REPORT-tier
exact numeric illustration of what the proven det-identity *means*.

**Matched-pair spectrum brick — positivity ⟹ real spectrum, machine-verified (cy420, 20th module).**
`forge/WitnessSpectrumReal.lean` (axiom-gate GREEN, `checked=8 clean=8 bad=0`) proves in Lean the punchline of
the matched pair, upgrading it from REPORT to HANDLE. For the two concrete `D'` matrices:
`det(D'_psd − s·I) = −s³ + 3s` (PSD witness `T ⪰ 0`) and `det(D'_ind − s·I) = −s³ − s` (indefinite witness).
Then `DpPSD_spectrum_real` exhibits the three real roots `{0, √3, −√3}` (PSD ⟹ **full real spectrum**), while
`DpInd_only_real_root_zero` proves the indefinite `D'` has `0` as its *only* real spectral point (`s²+1 > 0`),
its other eigenvalues `±i` being non-real (**indefinite ⟹ spectrum leaves the line**). So positivity of the
Weil form — CM's even-simple — is machine-verified to be exactly the condition that pins the spectral
realization's zeros onto the critical line. (The general abstract theorem "symmetric-in-a-PD-metric ⟹ real
spectrum" is a heavier deferred target — Mathlib `B^{1/2}`-conjugation + `charpoly_units_conj` + real-symmetric
splitting; pieces partially present.)

**Module inventory: 20** (19 + `WitnessSpectrumReal`).

**General abstract payoff — positivity ⟹ real spectrum (cy421–422, 21st–22nd modules).**
`forge/ConjSymmSelfAdjoint.lean` (axiom-gate GREEN, 3 theorems) proves the abstract, *dimension-general* form of
the CM payoff, sqrt-free (Mathlib has no matrix square root; positivity is instead expressed in Gram/factored
form `B = Pᵀ·P`, which is exactly how Weil positivity is witnessed):
- `conj_symm_of_selfadjoint`: if `B = Pᵀ·P` (`P` invertible) and `M` is B-self-adjoint (`Mᵀ·B = B·M`), then
  `P·M·P⁻¹` is symmetric;
- `charpoly_conj_eq`: `(P·M·P⁻¹).charpoly = M.charpoly` (similarity);
- `real_spectrum_of_selfadjoint_gram`: therefore `M.charpoly = ∏ᵢ (X − C eᵢ)` for real `eᵢ` — **M's spectrum
  is entirely real** (via Mathlib `Matrix.IsHermitian.charpoly_eq`).

This is the abstract "positivity ⟹ real spectrum" mechanism at full generality. **Honest scope (verify-before-
ship):** the master theorem's hypothesis `T·ξ = 0` makes the Weil form `T` *singular* (PSD, 1-dim kernel `⟨ξ⟩`),
not positive-*definite*; so `real_spectrum_of_selfadjoint_gram` (which needs `P` invertible ⟺ `T` PD) delivers
real spectrum for `D''` on the **quotient** `E/⟨ξ⟩`, where the induced Weil form is genuinely PD — *not* for `D'`
on the full space. It supplies the payoff half of "Weil positivity ⟹ RH"; wiring it to `D''` still needs the
deferred `D' → D''` quotient passage (`det0`, §E — Mathlib quotient-determinant infra). So the abstract chain is
now: master theorem (even-simple ⟹ `D''` self-adjoint in `T` + `ξ̂`-realization) + this payoff (T PD ⟹ real
spectrum) — closed **modulo** the quotient passage and the external even-simple positivity itself.

**Module inventory: 22** (+ `ConjSymmSelfAdjoint`).

**Positivity as a natural hypothesis — PD ⟹ Gram, and the capstone (cy423, 23rd module).**
`forge/PosDefGram.lean` (axiom-gate GREEN) removes the "already-factored" restriction:
- `posdef_gram_factorization`: `A.PosDef ⟹ ∃ P, IsUnit P.det ∧ A = Pᵀ·P` — built from Mathlib's
  continuous-functional-calculus square root `CFC.sqrt` (search-before-build found it in
  `Analysis/Matrix/Order.lean`; `P = √A` is Hermitian, `√A·√A = A`, invertible since `(det √A)² = det A ≠ 0`);
- `real_spectrum_of_selfadjoint_posdef` (**capstone**): `B.PosDef ∧ Mᵀ·B = B·M ⟹ M.charpoly = ∏ᵢ(X − C eᵢ)`,
  real `eᵢ` — self-adjoint in *any* positive-definite form ⟹ entirely real spectrum.

So the abstract "positivity ⟹ real spectrum" payoff is complete at its natural level of generality. **Honest
scope (unchanged):** applies to `D''` on the quotient `E/⟨ξ⟩` (induced Weil form PD), since the full-space `T`
is only PSD (`T·ξ = 0`). The remaining wire is the concrete **compression** brick (change of basis `Q = [ξ | W]`:
`Q⁻¹D'Q = [[0,b],[0,A]]`, `QᵀTQ = [[0,0],[0,T_W]]` with `T_W` PD, self-adjointness descending to `Aᵀ·T_W = T_W·A`
— then the capstone applies to the `(n−1)`-block `A`, and `det(D'−sI) = −s·det(A−sI)` is all-real). This
*sidesteps* the deferred quotient-determinant infra (only a change of basis, no quotient `LinearMap.det`).

**Module inventory: 23** (+ `PosDefGram`).

**Compression back-half — block-triangular ⟹ real spectrum (cy424, 24th module).**
`forge/RealSpectrumBlock.lean` (axiom-gate GREEN): `real_spectrum_blocktri` — for a block upper-triangular real
matrix `fromBlocks M₁₁ M₁₂ 0 A` with `A` self-adjoint in a PD form `T_W` and `M₁₁` already real-spectrum, the
whole `charpoly` splits into real factors; and `real_spectrum_blocktri_zero`, the CM specialization with
`M₁₁ = 0₁ₓ₁` (the `ξ`-kernel eigenvalue `0`). The compression mechanism was numerically verified end-to-end on
the cy419 PSD witness first (`Q⁻¹D'Q` block-triangular, `QᵀTQ` block-diagonal with PD `T_W`, `Aᵀ·T_W = T_W·A`,
`spec(A) = {±√3}`). So once `D'` is in `ξ`-basis block form, it has all-real spectrum. The remaining
**front-half** (fiddly, deferred): construct the change of basis `Q = [ξ | complement]`, prove `D'` equals the
block form there, and prove `T_W` PD via a *new* "PSD restricted off its radical is PD" lemma (Mathlib's
`PosDef.submatrix` is PD→PD, not PSD→PD). Still finite-N conditional on the external even-simple positivity and
the `N→∞` convergence — this is self-containment polish, not movement on the crux.

**Module inventory: 24** (+ `RealSpectrumBlock`).

**Front-half part-1 — PSD compressed off its radical is PD (cy425, 25th module).**
`forge/PosDefCompress.lean` (axiom-gate GREEN): `posDef_compress` — `T.PosSemidef` and `∀v, T·(W·v)=0 → v=0`
(the complement `W` meets `T`'s radical only at `0`) ⟹ `Wᴴ·T·W` is positive-**definite**. This is precisely the
PSD→PD "off the radical" fact that Mathlib's `PosDef.submatrix` (PD→PD only) lacks — the piece flagged in cy424.
Proof composes `conjTranspose_mul_mul_same` (`Wᴴ·T·W` PSD) with the PSD radical lemma
`PosSemidef.dotProduct_mulVec_zero_iff`. So all the *mathematical* content of the compression is now proven
(this + `real_spectrum_blocktri` + capstone + gram factorization); only the **change-of-basis assembly**
(construct `Q = [ξ | W]`, derive the two block identities) remains — pure bookkeeping. Still finite-N conditional
on the external even-simple positivity and `N→∞` convergence.

**Module inventory: 25** (+ `PosDefCompress`).

**Compression assembly connector (cy426, 26th module).**
`forge/CompressionAssembly.lean` (axiom-gate GREEN): `real_spectrum_of_similar_block` — if `Dp` on `Fin 1 ⊕ k`
is similar (`Q⁻¹·Dp·Q`) to the CM block form `fromBlocks 0 M₁₂ 0 A` with `A` self-adjoint in a PD `T_W`, then
`Dp.charpoly = ∏(X − C eᵢ)` (real spectrum). This wires the whole chain together (charpoly-similarity +
`real_spectrum_blocktri_zero`, PD block from `posDef_compress`). The **only** remaining gap is the standard
change-of-basis *existence* step — that the actual `D'` (`D'ξ=0`, self-adjoint in the PSD Weil form) is similar
to this block form via `Q = [ξ | complement]` (Mathlib `Basis.mkFinCons`). Pure basis bookkeeping, no remaining
mathematical content; still finite-N conditional on the external even-simple positivity and `N→∞` convergence.

**Module inventory: 26** (+ `CompressionAssembly`).

**Basis-descent brick + capstone (cy434, 27th module — closes the finite theorem).**
`forge/BasisDescent.lean` (axiom-gate GREEN, `checked=3 clean=3 bad=0`, `{propext, Classical.choice,
Quot.sound}`) — two theorems that **discharge the change-of-basis gap** left open above:

* `real_spectrum_of_radical_basis` — for **any** invertible `Q` whose first column is the radical vector `ξ`,
  it *derives* both the CM block form `Q⁻¹·Dp·Q = fromBlocks 0 M₁₂ 0 A` **and** the compressed self-adjointness
  `Aᵀ·T_W = T_W·A` from the natural hypotheses (`Dp·ξ=0`, `T` symmetric-PSD, `Dpᵀ·T=T·Dp`, `T·ξ=0`, `T_W`
  PosDef), then feeds `real_spectrum_of_similar_block`. The block form is automatic because `Dp·ξ=0` kills the
  whole `inl`-column of `B:=Q⁻¹DpQ`; the self-adjointness descends because `Dp` self-adjoint in `T` transports
  to `B` self-adjoint in the Gram form `QᵀTQ`, whose `inl` row/column vanish (`T·ξ=0`, `T` symmetric), leaving
  the `₂₂`-block identity `Aᵀ·T_W=T_W·A`. **This is exactly the "basis `{ξ, w_j}`, `D'` triangular" step** that
  was previously deferred — no more "pure basis bookkeeping" gap.
* `real_spectrum_of_weil_radical` (**capstone**) — the *faithful* finite form of **Connes–van Suijlekom
  Theorem 5.6(i)**: it assumes only the primitive CvS hypotheses — `T` PosSemidef, a **1-dimensional radical**
  (`hrad`: every `T`-null vector is a scalar multiple of `ξ`), `Dp·ξ=0`, and `Dp` self-adjoint in `T` — and
  concludes real spectrum. The positive-definiteness of the compression is **not assumed**: it is produced
  internally from `hrad` (+ `Q`-injectivity) via the `posDef_compress` brick, matching the paper's theorem
  statement rather than a weaker conditional. Magnitude-free (sign/kernel/block structure only; **no `ε_λ`
  dependence**). Still finite-`N`, conditional on the two external analytic inputs (even-simple positivity for
  the actual `QW_λ`, and `N→∞` Hurwitz convergence) — both untouched by this brick, both confirmed open by
  Connes 2026 §6.6.

**Module inventory: 27** (+ `BasisDescent`, 2 theorems). The abstract *"Weil positivity ⟹ real spectrum"*
reduction is now **complete end-to-end at finite `N`** — the finite core of the named CvS Thm 5.6(i) is fully
machine-verified, with the change-of-basis discharged. What remains is external analytic math, not formalization.
