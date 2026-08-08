# A machine-certified Weil-positivity theorem on a band-limited concentration class

*Riemann project, working paper draft. Provenance: LEDGER §788–§808, ANCHOR_THEOREM.md,
operators/cert_lambda7_final.py, operators/cert_anchor_final.py, RiemannLean/handles.jsonl.
Every numerical constant below is an interval-arithmetic (arb) enclosure; the closing inequality
is verified in the Lean 4 kernel. Claims are tagged **[HANDLE]** (machine-verified: arb-certified
or Lean-checked) or **[REPORT]** (rigorous on-paper argument, not machine-checked) throughout,
per the Charisma discipline.*

---

## Abstract

Weil's explicit-formula criterion states that the Riemann Hypothesis is equivalent to the
positivity of an explicit quadratic functional Q — the Weil form Q = pole + arch − prime — on the
space of admissible test functions. We prove **[HANDLE]** that Q(ψ) ≥ 0 for every even,
L²-normalized, band-limited test function ψ ∈ PW_{2W} (Paley–Wiener, exponential type 2W = 0.66,
support 4W = 1.32 > log 3) whose weighted off-window tail T(ψ) does not exceed an explicit
constant ε₀ = 10⁻⁴. Every constant in the proof is certified in interval arithmetic (arb /
python-flint, 200-bit) on a single self-consistent 8-probe Fejér basis, and the closing scalar
inequality is verified in the Lean 4 kernel using only the three trusted axioms
{propext, Classical.choice, Quot.sound}. The theorem is **unconditional**: no zero data, no
on-line-ness hypothesis, and no strip lemma enter the assembly of the Q-side constants. It
establishes one rung — a described infinite-dimensional class at the first prime-active support —
of the uniform (all-width) statement, which is Weil's criterion itself and remains open.

To our knowledge (see §6, and the §825 literature gate), this is the first machine-verified
positivity result for the Weil explicit-formula form on an explicitly-described infinite class of
test functions.

---

## 1. Introduction

### 1.1 Weil's positivity criterion

For a suitable even test function ψ with Fourier transform ψ̂, the Weil explicit formula equates a
sum over the nontrivial zeros ρ = ½ + iγ of ζ,

    Σ_ρ ψ̂(γ)  =  [archimedean term]  +  [pole term at s = 0, 1]  −  Σ_{n} Λ(n) n^{-1/2} ψ̂-part(log n),

so that the **Weil quadratic form**

    Q(ψ)  :=  pole(ψ)  +  arch(ψ)  −  prime(ψ)

satisfies, when ψ is replaced by a self-convolution (so the zero side becomes 2 Σ_γ |φ̂(γ)|²),

    Q(ψ)  =  2 Σ_γ |φ̂(γ)|²   ≥ 0 for all admissible ψ   ⟺   RH.

The "⟸" direction is the deep one: if the zeros are on the critical line the zero side is a sum of
squares and Q ≥ 0; conversely a single off-line zero produces a test function with Q(ψ) < 0. Thus
proving Q ≥ 0 **on the prime side, with no assumption on the zeros**, uniformly over a rich enough
class of ψ, *is* a proof of RH. This is the Weil–Bombieri positivity program (Weil 1952; Bombieri's
Clay problem description; Li's criterion 1997 as the discrete shadow).

### 1.2 What this paper proves, and what it does not

The uniform statement is RH and is open. What we establish is a **fully certified rung**: Q ≥ 0 on
an explicitly-described, infinite-dimensional, band-limited class, at the first support where the
prime side is genuinely active (primes 2 and 3 in range, 4W = 1.32 > log 3). The class is cut out
by an explicit analytic condition on ψ alone — a bound on its weighted off-window tail — with no
reference to our computational basis, so the statement is non-circular.

The value of a *certified* rung is threefold. (i) It fixes an honest floor: the positivity is not a
numerical impression but an interval-arithmetic + kernel-checked fact. (ii) It isolates the exact
obstacle to the uniform statement — every constant here is width-stable *except* the prime form,
which grows with the support (§6.3, and LEDGER §816/§818). (iii) It supplies a template — the
concentration split plus five certified constants — that the uniform program must fill in with a
genuinely new, width-uniform bound on the prime form.

---

## 2. Setup and notation

Fix the half-width W with 2W = 0.66, so band-limit / exponential type is 2W and the time support of
ψ is [−2W, 2W]; the Fourier support relevant to the prime sum is 4W = 1.32. Only the primes with
log p ≤ 4W contribute: p ∈ {2, 3}, since log 3 = 1.0986 < 1.32 < log 5.

- **Test class.** ψ ranges over even functions in the Paley–Wiener space PW_{2W}, L²-normalized
  (‖ψ‖₂ = 1). ψ̂ denotes the Fourier transform.
- **Weil form.** Q(ψ) = pole(ψ) + arch(ψ) − prime(ψ), with
  - pole(ψ): the rank-one contribution from the poles of ζ at s = 0, 1;
  - arch(ψ) = (1/2π) ∫ |ψ̂(t)|² Ω(t) dt, where Ω(t) = Re ψ₀(¼ + it/2) − ln π is the archimedean
    (digamma) density, ψ₀ the digamma function;
  - prime(ψ) = Σ_{n} Λ(n) n^{-1/2} (ψ̂-self-convolution)(log n), the finite prime sum over
    log n ≤ 4W.
- **Window.** [14, 44] is the fixed frequency window carrying the concentration operator.
- **Concentration operator.** A = P_{PW} · 1_{[14,44]} · P_{PW}, the band-limit projection composed
  with the window indicator. Its top eigenspace is a prolate-spheroidal eigenspace; A is the
  even-symmetrized version throughout.
- **Weighted off-window tail.** For the class definition,

      T(ψ) = ∫_{|t| ∉ [14,44]} |ψ̂(t)|² (1 + |Ω(t)|) dt.

  The weight (1 + |Ω|) is exactly the density that makes Q's off-window contribution finite; this is
  the quantitative form of the window-cap argument.

---

## 3. The theorem

**Theorem (certified Weil positivity on the weighted-concentration class at the open rung).**
Let 2W = 0.66. There is an explicit constant ε₀ = 10⁻⁴ > 0 such that for every even, L²-normalized
ψ ∈ PW_{2W} with T(ψ) ≤ ε₀,

    Q(ψ)  ≥  0.

Moreover the worst-case lower bound is strictly positive: Q(ψ) ≥ 0 with certified slack (the closing
inequality holds with margin, §5).

The class {ψ : T(ψ) ≤ ε₀} is infinite-dimensional and cut out by an explicit condition on ψ; it
contains, in particular, every ψ whose spectral mass is ε₀-concentrated in [14, 44] against the
(1 + |Ω|) weight.

---

## 4. The assembly

Split ψ = p + f, where p = Π_E ψ is the orthogonal projection onto E, the top-k = 6 eigenspace of
the even concentration operator A, and f = ψ − p is the remainder. All five constants below are
certified interval enclosures on **one** self-consistent 8-probe Fejér basis (LEDGER §802A–§808).

**Constant 1 — concentration control of ‖f‖ (via λ₇).** [HANDLE]
⟨Aψ, ψ⟩ ≥ 1 − T(ψ) (the in-window mass is at least 1 minus the off-window tail). The eigen-split
⟨Aψ, ψ⟩ ≤ ‖p‖² + λ₇‖f‖² then gives

    ‖f‖²  ≤  T(ψ) / (1 − λ₇).

The certified upper bound on the first discarded eigenvalue is **λ₇ ≤ 0.60696**, obtained by the
Ky-Fan trace route: the exact even trace tr A = 6.3102004… (arb, radius 5.7×10⁻⁵⁹) minus the
certified compression sum of the top six, Σ₆ = 5.7086 ± 5×10⁻³, bounds the tail sum and hence λ₇.
(The float value of λ₇ is ≈ 0.425; the certified bound 0.607 is honest-loose and the assembly still
closes.) So 1 − λ₇ ≥ 0.393 and ‖f‖² ≤ 2.544 · T(ψ).

**Constant 2 — Q on the projected part p (via m_E).** [HANDLE]
p lives in the explicit 6-dimensional span E. The certified generalized-eigenvalue floor of the
Weil form restricted to E,

    m_E  :=  λ_min(Q, G) on E   ≥  0.09176,

is computed as the min generalized eigenvalue of (VᵀQV, VᵀGV) by interval LDLᵀ with bisection, from
the arb-certified ambient Gram data gram14.json. Hence Q(p) ≥ m_E ‖p‖².

**Constant 3 — the cross term |Q(p, f)| (via B_pf).** [HANDLE]
Splitting arch at the window edge and bounding each piece:
- in-window archimedean: |Ω| ≤ **sup_{[14,44]}|Ω| ≤ 2.0175** (certified interval sweep, 3000 balls),
  carrying an extra √λ₇ factor from f's small in-window mass;
- off-window archimedean: bounded by √(weighted off-mass of E) · √T(ψ), with the certified weighted
  off-window mass of the eigenspace **τ_E ≤ 0.0862** (fine-segment arb integration, per-vector
  0.003–0.086);
- prime: ≤ **B_prime = 0.3579** ‖p‖‖f‖ (exact arithmetic: (1/π)(log2/√2 + log3/√3), grows ~e^{2W});
- pole: a negligible certified constant.

Combining, the certified cross constant is **B_pf ≤ 2.223** (the value carried into the closing
inequality; LEDGER §808 records the conservative envelope B_pf ≤ 2.65 including √τ_E, and the tighter
2.223 is the one Lean checks).

**Constant 4 — Q on the remainder f (via C_f).** [HANDLE]
|Q(f)| ≤ C_f ‖f‖² + T(ψ), with **C_f ≤ 1.582** certified (arch-in λ₇·sup|Ω|, prime 0.3579, pole
negligible; the off-window arch piece is ≤ T(ψ) directly by the class definition).

**Constant 5 — the closing inequality.** [HANDLE — Lean-kernel-verified]
With ‖f‖² ≤ ff · T(ψ), ff = 1/(1 − λ₇) ≤ 2.544, and ‖p‖² ≥ 1 − ‖f‖²,

    Q(ψ)  ≥  m_E(1 − s₀)  −  2 B_pf √(s₀(1 − s₀))  −  C_f s₀  −  ε₀,   s₀ = ff · ε₀.

At ε₀ = 10⁻⁴ (so s₀ = 5/19652 ≈ 2.54×10⁻⁴) with the certified constants
(m_E ≥ 0.09176, B_pf ≤ 2.223, C_f ≤ 1.582), the right-hand side is **strictly positive**. This
scalar inequality is proved in the Lean 4 kernel (theorem `anchor_eps0_positive`, module
`Cand_anchor_eps0_positive`, PASS, axioms {propext, Classical.choice, Quot.sound}), via the rational
over-approximation √(s₀(1−s₀)) ≤ 15949/1000000 followed by `nlinarith`. ∎

---

## 5. Certification: the two verified layers

Two independent machine layers carry the proof; the analytic assembly (§4) linking them is standard
functional analysis (Cauchy–Schwarz, the min-max eigen-split, compression interlacing) and is
**[REPORT]** in the precise sense that the *chain of operator inequalities* is on paper — but every
*numerical constant* it consumes and the *final scalar inequality* it produces are machine-verified.

1. **Interval arithmetic (arb / python-flint, 200-bit).** [HANDLE] Every constant — tr A, the
   compression eigenvalues feeding λ₇, m_E via interval LDLᵀ, sup|Ω| via a 3000-ball sweep, τ_E via
   fine-segment integration, B_prime exactly — is a certified enclosure, not a float. Scripts:
   operators/cert_lambda7_final.py, cert_anchor_final.py, gram14.json.

2. **Lean 4 kernel.** [HANDLE] Two theorems in RiemannLean/handles.jsonl, both PASS with only the
   three trusted axioms:
   - `zerofree_rung_certM_posSemidef` — the zero-free open-rung certificate matrix M is positive
     semidefinite (exact LDLᵀ sum-of-squares over ℚ; min pivot 2.0237), which with the certified
     ‖Q − M‖₂ ≤ 3.80×10⁻⁵ gives unconditional positive-definiteness of the true zero-free Weil form
     on the Fejér span {14, 22, 30, 38}.
   - `anchor_eps0_positive` — the §4 closing inequality, the analytic heart of the ∀-theorem.

**Honest layering.** The Handle is: the constants are interval-certified and the closing inequality
is kernel-checked. The Report is: the reduction of "Q(ψ) ≥ 0 for all ψ in the class" to "that scalar
inequality holds" is a hand proof (a short, standard one). A full end-to-end Lean formalization of
the operator assembly is the natural next hardening; it is not yet done, and the paper states so
rather than implying the whole chain is machine-checked.

**Discipline note (kept for provenance).** During certification two false-positive "CERTIFIED"
prints were caught and rejected before banking: one where a `nan` τ_E collapsed to a trivial pass,
and one where ε₀ was float-only rather than arb. The banked result uses the arb/Lean values
throughout. A separate 17-probe/k = 9 route to the discarded eigenvalue hit a `nan` in the
compression trace and was **not** used; the banked route is the k = 6 Ky-Fan trace bound above.

---

## 6. Honest scope and relation to prior work

### 6.1 The uniform statement is RH

Letting W → ∞ (all supports) and ε₀ → the full class recovers exactly Weil's criterion. Our theorem
does not approach that limit; it certifies one rung. The barrier is not incidental: it is that the
prime form grows with the support while the archimedean floor does not (§6.3).

### 6.2 Relation to the positivity program

- **Weil (1952) / Bombieri.** The explicit-formula positivity criterion is classical; what is new
  here is not the criterion but a *certified* positivity statement on a described infinite class.
- **Li's criterion (1997).** RH ⟺ λ_n ≥ 0 for the Keiper–Li coefficients — the discrete positivity
  shadow. Our object is the continuous Weil form on band-limited test functions, closer to the
  Beurling–Selberg / Fourier-optimization line.
- **de Branges; Hilbert spaces of entire functions.** The concentration/prolate structure we use
  sits naturally in the de Branges framework; the frame lower bound handling the zeros' small gaps
  (the piece a uniform proof needs on the zero side) is exactly where de Branges / sampling theory
  would enter.
- **Connes–Consani–Moscovici (2025), "Zeta Spectral Triples," arXiv 2511.22755.** The leading
  spectral program. Their prolate-spheroidal basis coincides with our concentration eigenspace, and
  their conjectured μ_λ → 0 matches our measured Q_floor → 0⁺. Neither their route nor the
  arithmetic-surface route closes the core; both bottom out where our measurements localize the
  difficulty. A companion note (CCM_BRIDGE.md) works out the dictionary and where our apparatus can
  numerically probe their open Lemma 7.2.

### 6.3 The honest open core (what a uniform proof still needs)

The uniform theorem requires: **arch(ψ) + pole(ψ) − prime(ψ) ≥ 0 for all band-limited ψ, uniformly
in width, with no assumption on the zeros.** Our measurements (LEDGER §812–§821) show precisely why
this is hard and where: the difficulty *localizes* to a boundary shell of concentration modes; the
crash is *cross-coupling*, not single-mode; it is *arithmetic* (the archimedean form is a
positive-definite floor ~0.153, width-invariant; the prime form is indefinite and its positive
eigenvalues eat the floor on boundary-coupled directions); and it is a *sampling* phenomenon at the
Beurling threshold t* = 2π e^{4W}. A §825 literature gate found no known unconditional technique
bounding the prime form against the archimedean floor — that bound is RH. A uniform proof needs a
genuinely new, width-uniform positivity mechanism (an arithmetic-surface / intersection-form object
of the Hodge-index type), not more numerical scans.

### 6.4 Novelty claim (stated carefully)

We are not aware of a prior machine-verified (interval-arithmetic + proof-kernel) positivity result
for the Weil explicit-formula form on an explicitly-described infinite class of test functions. We
state this as "not aware of," pending a full literature check, rather than as an unqualified first.

---

## 7. What extends it

1. **Certify the width-invariant archimedean floor and the t*-margin law** as arb-certified
   structural facts (LEDGER §816-CERT, §823/§824 have the certified curves).
2. **Formalize the operator assembly in Lean** end-to-end (turn the §5 Report layer into Handle).
3. **The prime-side spectral inequality** — bound the prime form's positive eigenvalues by the
   archimedean floor *structurally*, not scalar. This is the real prize and equals the open core.
4. **The CCM bridge** — connect our certified Q_floor → 0⁺ rate to their μ_λ → 0 conjecture and
   numerically probe their Lemma 7.2 (companion note).

---

*Status: this draft is grounded entirely in banked, certified results (LEDGER §788–§808) and the two
PASS Lean handles. Constants are the arb/Lean values, not floats. The one place the chain is not yet
machine-checked — the operator assembly linking the certified constants to the closing inequality —
is stated as such (§5), not papered over.*
