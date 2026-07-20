# Vector D — function-field / Hodge-index transport: scoping (2026-07-19)

**Context (search-before-build).** The analytic toolkit is fully mapped in `OBSTRUCTION_MAP.md`: 7 routes,
all reduce to ★ (A+B ≥ P), all shown to be a *global* Diophantine constraint. This session's compute phase
(LEDGER §10–§17) independently re-derived routes 4 (A+B indefinite, low-n negative modes) and 5 (every
operator certificate reduces to Weil). Vector D is the ONE route NOT in that map because it is **geometric**,
not analytic. This doc scopes it honestly: it is the deepest frontier in the subject (Deninger; Connes–
Consani F₁), and no overnight computation cracks it. The value here is a precise dictionary + the exact gap.

---

## 1. Why the function-field case is PROVEN (the mechanism to transport)

For a smooth projective curve C over F_q, the Hasse–Weil zeta zeros satisfy |α_i| = √q ("RH for C"). Proof
(Weil 1948; Grothendieck/Deligne cohomological): work on the **surface S = C ×_{F_q} C**.
- The zeros α_i are eigenvalues of **Frobenius F acting on H¹(C)** (ℓ-adic or crystalline cohomology).
- On the surface S, divisor classes have an **intersection pairing** ⟨·,·⟩. The **Hodge Index Theorem**:
  the intersection form on the Néron–Severi group NS(S) ⊗ ℝ has signature **(1, ρ−1)** — exactly one
  positive direction (the hyperplane/ample class), all others negative.
- Restricting to the **primitive** part (⊥ to the trivial correspondences: the two fiber classes and the
  diagonal), the form is **negative definite**. Applied to the graph of Frobenius Γ_F, the
  **Castelnuovo–Severi inequality** (a positivity from that definite form) gives |trace F|² ≤ (deg terms),
  which forces |α_i| = √q. **Positivity of a definite intersection form ⟹ RH.**

The engine is: an **arithmetic surface** + an **intersection form** + a **Hodge/Lefschetz sign theorem**
that makes a specific quadratic form definite. That definiteness IS the RH positivity.

---

## 2. The number-field dictionary (what each object should become for ζ(s))

| function field (curve C / F_q)          | number field (ℚ, ζ(s)) — the sought analog                        |
|-----------------------------------------|-------------------------------------------------------------------|
| curve C over F_q                        | Spec ℤ "compactified" ( + the archimedean place ∞ ) = "Spec ℤ̄"    |
| base field F_q                          | the "field with one element" **F₁** (does not classically exist)  |
| surface S = C ×_{F_q} C                 | **Spec ℤ ×_{F₁} Spec ℤ** — the missing object                      |
| Frobenius F on H¹(C)                     | a self-adjoint/scaling operator with spectrum = {γ : ζ(½+iγ)=0}   |
|                                          | (Hilbert–Pólya; Connes' scaling action on the adele class space)  |
| intersection pairing on NS(S)           | **our τ = A + B − P** (the Weil explicit-formula quadratic form)   |
| Hodge Index / Castelnuovo positivity    | τ ⪰ 0 (= RH) — the definiteness we keep certifying numerically     |
| deg = q on the diagonal                  | the archimedean/pole normalization (B block)                       |

**Our τ is already sitting in the right slot:** it is the arithmetic analog of the intersection form whose
definiteness would prove RH. This session's structural map of τ is a map of that intersection form.

---

## 3. The exact missing ingredient (the gap, stated precisely)

The transport fails at ONE place, and it is the central open problem of the whole geometric program:

> **There is no base for Spec ℤ.** A curve C lives *over* F_q, which lets you form the self-product surface
> S = C ×_{F_q} C carrying the intersection theory that powers the proof. Spec ℤ has no base field, so
> **Spec ℤ ×_{F₁} Spec ℤ and its intersection/Hodge theory do not (yet) exist.** Everything else in the
> dictionary has a candidate; this one object does not.

Two active programs try to supply it (both referenced so this integrates, not duplicates):
- **Connes–Consani F₁ geometry** — the *arithmetic site* and *scaling site*; a Riemann–Roch on the adele
  class space; our τ comes from the Connes–van Suijlekom (`cvs`) construction in exactly this lineage
  (see [[connes-decisive-diagnostic]], `operators/connes_decisive_v3.py`). They obtain an H¹-like space and
  a trace formula, but **the intersection/Hodge positivity — the analog of the definite form — is the piece
  still missing.** That missing piece is literally ★.
- **Deninger's cohomology** — a conjectural infinite-dimensional cohomology + a flow whose "Lefschetz trace"
  is the explicit formula; the sought positivity would come from a polarization (Hodge–Riemann). Also open.

**Consequence:** Vector D does NOT bypass ★ — it would *derive* ★ from a geometric definiteness theorem
instead of proving it analytically. That is why it is attractive (it is the route whose positivity is also
global/geometric, matching RH's character per OBSTRUCTION_MAP's "common root") and why it is hard (it needs
a new geometric object, not a new estimate).

---

## 4. What is actually approachable (honest triage)

- **NOT approachable overnight / by us alone:** constructing Spec ℤ ×_{F₁} Spec ℤ or Deninger's cohomology.
  This is frontier research, decades deep, and a script cannot produce it.
- **Approachable — concrete probes that could reveal whether the geometric backbone is present in our τ:**
  1. **Hard-Lefschetz / sl₂ test:** the Hodge Index Theorem is powered by an sl₂ (Lefschetz L, Λ, grading H)
     acting on cohomology, whose primitive decomposition + Hodge–Riemann relations give the signature. TEST:
     does our τ (or A+B) carry a compatible sl₂ — a raising operator on the frequency grading under which
     the spectrum organizes into Lefschetz strings? If yes → a real geometric echo to exploit. If no → the
     naive geometric structure is absent and Vector D needs the full F₁ object (informative either way).
  2. **Frobenius-unitarity shadow:** in the function-field case positivity ⟺ Frobenius is unitary after
     scaling (eigenvalues on |α|=√q). Our zeros-operator analog (the Connes scaling operator, not τ) should
     be self-adjoint; τ's definiteness is its form-level shadow. Check the relationship numerically.
  3. **Inertia vs Hodge signature:** compare τ's / A+B's inertia growth (n_neg ~ 1.8 log c, §12) to any
     Hodge-index prediction (1, ρ−1) for a candidate arithmetic-surface Picard number ρ. (Prior read: no
     clean match — A+B is mostly-positive, Hodge form is mostly-negative — but worth stating precisely.)

**Recommended first probe:** #1 (hard-Lefschetz/sl₂), because it directly asks "is the algebraic engine of
the function-field proof present in our operator?" — a crisp yes/no with either outcome informative.

---

## 5. Honest status
Vector D is the correct *character* of route (global/geometric positivity) but requires a geometric object
that is the subject's central open problem. We cannot build that object. We CAN test whether our τ already
carries its algebraic shadow (probe #1). If it does, that is a genuine new handle; if it does not, we have
ruled out the naive geometric transport and confirmed the F₁ object is genuinely required. Everything stays
REPORT tier. This is scoping, not progress on the proof.

---

## 6. Probe #1 result (2026-07-19): hard-Lefschetz backbone is ABSENT in τ
Ran `operators/lefschetz_probe_D1.py` (results in-line). Both tests NEGATIVE at c=14,50,200:
- **Hodge–Riemann sign pattern:** A+B's negative modes have MIXED frequency parity (odd/even =
  2/2, 4/3, 6/4) — not the single-parity pattern an sl₂ weight structure would force.
- **Lefschetz self-adjointness:** the +1 frequency shift is NOT τ-self-adjoint (‖τL−Lᵀτ‖/‖τL‖ ≈ 1.38).

**Conclusion:** the naive hard-Lefschetz/sl₂ algebraic engine of the Hodge Index Theorem is NOT already
present in τ. Vector D cannot shortcut by reading a geometric structure off the operator; it genuinely
requires constructing the missing F₁ object (Spec ℤ ×_{F₁} Spec ℤ + its intersection theory). This confirms
the gap in §3 is real and not an artifact of not-having-looked. Informative negative; Vector D remains the
correct *character* of route but is frontier-hard, not overnight-tractable.
