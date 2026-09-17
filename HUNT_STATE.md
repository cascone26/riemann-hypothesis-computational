# Riemann Hypothesis Hunt — State Summary (2026-08-10)

**Compiled from:** OBSTRUCTION_MAP, BARRIER_MAP, RH_DIFFICULTY_CONSERVATION_MAP, RH_WALL_ANATOMY, RH_CORE_LOCALIZATION, MASTER_SYNTHESIS, DENINGER_RELATED_WORK, LEDGER §884–§981, all verified/double-checked.

**Purpose:** Sharp anti-backtrack reference. One line per result; no editorial. Read this before starting a session to avoid re-deriving or re-attempting dead ends.

---

## 1. ESTABLISHED (verified, don't re-derive)

- **Li coefficients λ_n > 0 for n = 1..3732 unconditional** (Theorem 3B, PROOF_SYNTHESIS); n ≤ 85 analytically proven, n = 86..3732 via pair-margin method on 100k LMFDB zeros.
- **de Bruijn–Newman Λ ≤ 0.176168** (certified, §499 MASTER_SYNTHESIS; beat published record 0.22; floor is real at Polymath15 mollifier exponent y₀).
- **Weil form τ = A + B − P ⪰ 0 ⟺ RH exactly** (five independent reformulations verified to same a=1 razor-zero margin; §470–474 MASTER_SYNTHESIS).
- **M_σ = arch − prime is Pontryagin Π₁: #neg(M_σ) = 1 always** (archimedean vs prime balance, ζ's simple pole absorbs the negative square, κ=1 cross-Selberg verified; §946–979 LEDGER).
- **Function-field RH proven via Hodge-index on C×C** (Weil 1948; gives positivity a priori from geometry; number-field analogue missing).
- **Archimedean density is Cauchy–Stieltjes, completely monotone** (dΩ/d(r²) = unconditional series from ψ₀').
- **Razor is off-diagonal / near-null coupling** (~2–3 dimensional prolate-plunge subspace, irreducibly degenerate, ~0.53 overlap; margin ≈ 1e-15, sub-interval-radius).
- **Krein string: phases arg ζ′(ρ) encode Möbius & Dirichlet-inverse coefficients universally** (Möbius squarefree-detection r=0.936, phase-only r=0.989; proven degree-1 ζ/χ and degree-2 Ramanujan Δ; §511–979 LEDGER).
- **κ = #pole-locations at s=1 (∈{0,1}), not order** (ζ,ζ²,ζ³ all κ=1; entire L-functions κ=0; degree-2 with one pole κ=1; robust at saturation).
- **Geometric containment: arch−prime PSD on pole^⊥ ⟺ RH** (the one negative square sits entirely in span(pole); §975–977 LEDGER, faithful to perturbation).

---

## 2. RULED OUT / DEAD ENDS (never re-attempt — reason each died)

| Approach | Death | Category |
|---|---|---|
| **Moment/Hankel (Li λ_n as Stieltjes moments)** | Hankel matrix fails PD at N=2 (det = −0.0037); (λ_n) is log-concave, opposite of moment sequence | K2: premise-false |
| **Beurling-Selberg pointwise majorant** | Prime comb Λ(k)/√k exceeds archimedean e^{-\|x\|/2}/(1−e^{−2\|x\|}) at every prime ≥5; ratio → ∞ | K2: pointwise fails |
| **Large sieve (global bounds on P)** | Smooth prime symbol on line gives zero-density (off-line only); cannot force zeros ON the line | Hits wall |
| **de Branges naive (H(E) reproducing kernel)** | A+B has 4 negative eigenvalues (§4 in OBSTRUCTION_MAP); only τ=A+B−P is PSD | K1: witness Conrey–Li |
| **Operator spectral certificates for τ⪰0** | All finite certificates (Gershgorin, SOS, tridiagonal) reduce to Weil positivity; τ⪰0 ≡ RH | K3: circular |
| **Scalar/operator floor−pnorm bound at corner** | Misalignment correction closes only ~75%; collapse nonlinear at knife-edge; §914/§919 CORE_LOCALIZATION | Fails locally |
| **Concentration via Szegő (band-limited)** | Diagonal symbol wrong-sign at corner; off-diagonal effect dominates; robust only away from edge | Robust-region-only |
| **Low-prime / first-few primes** | Margin is global cancellation over FULL prime set; cutting primes destroys the balance; §917/§925 CORE_LOCALIZATION | Global-only |
| **Two-term arch ≥ p=2 dominance** | Fails at 4W≥6; primes cannot be peeled; every rung requires all prior primes (§485 MASTER_SYNTHESIS) | Peeling fails |
| **Density-preserving (PNT smooth + error terms)** | Smooth part vacuous at window; entire content IS the error term ≡ RH; §926 CORE_LOCALIZATION | Circularity |
| **de Branges global (ω→0 in canonical system)** | Conrey–Li 2000 counterexample off-line zeros at 34th zero (W1 Davenport–Heilbronn satisfies axioms) | K1: witness realized |
| **Function-field / Hodge translation to ℤ** | Deninger 2022: no real-coefficient arithmetic Weil cohomology on number fields (parity obstruction: odd multiplicity symplectic vs even complexified) | K4: theorem-required |
| **Arb-certification of the razor margin** | Margin ~1e-15, below achievable interval radius (~1e-13 at best); §894 CORE_LOCALIZATION | Tool floor |
| **Prime-induction / finite-core reduction** | Knife-edge margin at p=2 (~1e-4) defeats every effective per-prime bound; collapse not self-similar (6.6e69× variance §14) | K4: ceiling |
| **Jensen/Laguerre-Pólya/total-positivity ladder** | TP2/hyperbolicity hold in smooth regime (wide functions, all heights) but fail at corner (dense regime where RH lives); strengthening does not gain slack | K4: partial axis |
| **Hilbert–Pólya (generic self-adjoint)** | Reality free, correspondence (spectrum=zeros) the hard half; conservation law bars single framework from having both halves free | Framework limit |
| **Trace formula / eigenvalue-sum identities** | Re-clothes explicit formula; positivity-cancellation ≡ RH; K5 identity; carries no new information | K5: identity |
| **Averaging / family-mean arguments** | Large-sieve forcing EVERY χ would prove GRH (a stricter statement); outlier-perturbations not bounded by variance; §62/§B10 BARRIER_MAP | K4: stricter |
| **DBN heat flow (λ≠0 deformation)** | Destroys Euler product; all RH-violators drop either EP or FE; deformation enters witness territory; §47d-e in-repo | K1: witness class |
| **Extremal / Viazovska-style magic function** | No sharp single extremal exists; near-null is irreducibly 23-dim degenerate subspace (grows with c); dual inherits primal degeneracy | K2: premise-false |
| **Zero-statistics (GUE pair correlation, Odlyzko moments)** | Montgomery pair-corr proven ASSUMING RH (circular); Alternative Hypothesis on half-integers consistent with all known data (Goldston et al.); GUE universal across unrelated families | K1+K3+K4: triple |
| **Approximate FE instead of exact Riemann** | Permits Turán's counterexample (partial-sum zeros); exact Γ-factor & growth regularity load-bearing; no middle ground | K2: weakening fails |
| **Naive m-function / de Branges E-function for κ-count** | E(z)=W(screw) has no UHP zero (winding 0); naive Nevanlinna-Gram sampling misses indefinite-metric structure; route open but not built | Not yet built |

---

## 3. GENUINELY OPEN / UNBUILT (not attempted, or attempted inconclusively)

- **Arithmetic-surface Hodge-index theorem (Spec ℤ ×_{F₁} Spec ℤ)**: Missing geometric positivity source; C–C have substrate candidates (2606.06604) but NO intersection-pairing / Hodge-index construction published; positivity crux unsolved. (The concrete keystone per §264–267 DIFFICULTY_CONSERVATION_MAP & §329 CORE_LOCALIZATION.)
- **Condensed-math escape (Clausen–Scholze solid/liquid real coefficients)**: Deninger's no-go is ordinary-real-specific; whether condensed/analytic real evades it is unexplored (no published work). Needs deep feasibility study.
- **Weyl m-function / Nevanlinna realization of M_σ**: Building the ζ-related m-function whose pole-count (argument principle) gives κ=1 directly; canonical-system duality open. The §974 "gap" in CORE_LOCALIZATION.
- **Full Π₁ canonical-system construction**: The de Branges system whose m-function inherits ζ's pole as a negative square and realizes κ-machinery; not built.
- **Higher RH-verification height (T > 3×10¹²)**: Platt–Trudgian bottleneck; pushing Λ below 0.176 requires new verified zero runs, large compute. Well-defined, incremental.
- **Better DBN mollifier barrier**: Polymath15's own open problem; {2,3,5,7,11,13} mollifier at ≤0.176 floor; next lever unexplored.
- **Full universality law for κ across degrees**: Started Selberg class (degree 1 & 2 cross-confirmed); degree-d generalization not complete. Candidate theorem.
- **Fractional integration / Riemann-Liouville frameworks**: Not systematically explored in this project.
- **Chromatic lens / algebraic-geometry re-encoding**: Eisenbud-Evans or similar structural reframings not attempted.

---

## Key Invariant (The Wall)

Every framework for RH gets **exactly one of {reality, correspondence} for free; the other is the hard half:**
- *Reality* = zeros on critical line (spectrum real / operator self-adjoint).
- *Correspondence* = zeros are ζ's zeros (built from primes / arithmetic).

Single frameworks cannot cross this divide (§289 DIFFICULTY_CONSERVATION_MAP). **The one exception:** function fields, where **geometry (Hodge-index)** supplies both halves at once via the surface C×C. Number fields lack that surface — RH lives in the gap.

---

**Next Session:** Start by reading this, then read CORE_LOCALIZATION §7 (the frontier ranking) and STATUS.md (project momentum). Do NOT attempt anything in section 2 without explicit justification of why-it-differs-this-time.

**ETA to frontier:** The missing keystone (Hodge-index / condensed-math / m-function realization) is a sustained research program, not accelerable compute. No session-level acceleration possible; long-game is the only game.

---
## FRESH FRONTS — 2026-08-10 redirect (riddle-mode; NOT yet attempted here)
Context: Jacob redirected from analytic grinding to riddle/reframe hunting. Anti-repeat spine confirmed the deepest
crux is the EULER PRODUCT (BARRIER_MAP: EP∧FE is the only property-pair with an empty RH-violator column; every
violator, e.g. Davenport–Heilbronn, drops exactly one). Two fresh fronts, neither in the 23 dead-ends:

- **FRONT A — ELEMENTARY / EULER-PRODUCT (Robin/Nicolas/Lagarias).** Analysis-free equivalents (σ(n)/n bounds,
  primorials, colossally-abundant numbers). Pure multiplicativity = the identified crux, naked. NEVER worked here
  (all 23 dead-ends are analytic/zeros-side). Committed first front. Survey of its frontier launched 2026-08-10.
- **FRONT B — QUASICRYSTAL / LEE–YANG (Dyson 2009; Kurasov–Sarnak 2020; Shaughnessy 2410.03673).** Zeta zeros as
  a 1D Fourier quasicrystal; build from stable/real-rooted polynomials. Connects to Tony's aperiodic-order work.
  CAUTION: this is the real-rootedness thread; its cousin (Jensen/Laguerre–Pólya ladder) is ALREADY ruled out
  ("TP2 fails at corner") — may be the same wall in Fourier clothing. Live #2, revisit after Front A.

Tony's current work (checked 2026-08-10): systems/epistemology + aperiodic order; DRY on direct RH math, but the
aperiodic-order structure is the bridge to Front B.

### Front A probe log (2026-08-10, hands-on, before frontier survey landed)
- Corrected intuition: primorials are NOT Robin-extremal (G→(6/π²)e^γ≈1.083, drifts AWAY). The extremal seq is
  COLOSSALLY ABUNDANT numbers (rigid exponents a_p decreasing in p, set by ε); there G→e^γ from below.
- DECOMPOSED the Robin margin at CA numbers: log(e^γ/G) = [combinatorial CUSHION −logC, rigid, +] +
  [residual]. Residual = −e_M + (loglogn−logP gap), where e_M = Mertens product error logM−(γ+log logP).
- FINDING (verify-before-ship; caught a loglog bug mid-probe): e_M is COMPARABLE to the cushion (both ~0.003–0.14),
  hence LOAD-BEARING. So the naive CA-decomposition does NOT bypass the zeros — the RH content localizes to the
  Mertens product error e_M, which is zero-governed (explicit formula for ψ(x)). Cushion vs |residual| crosses at
  exactly n≈5040 (Robin's threshold falls out). NOT a kill of Front A, but: the *naive* combinatorial decomposition
  reconnects to zeros via e_M. Open: does a cleverer combinatorial handle (Nicolas primorial form? extremal-number
  rigidity used differently?) avoid e_M? — deferred to the elementary-frontier survey (running).

### Front A — CLOSED as a bypass (2026-08-10, hands-on + cited survey)
DEAD. The elementary equivalents (Robin, Nicolas, Lagarias, Caveney–Nicolas–Sondow) do NOT bypass the zeros:
1. Hands-on: Robin margin = [rigid combinatorial cushion] − [Mertens error e_M]; Nicolas = pure e_M. RH-content = e_M
   = prime-counting error π(x)−li(x) = the zeros. (Cushion vs |residual| crosses at exactly n≈5040, Robin's threshold.)
2. Survey (cited): the natural combinatorial handle — a MONOTONICITY proof of Nicolas — PROVABLY contradicts Cramér's
   prime-gap conjecture (arXiv:1012.3613). Prime-gap oscillation kills monotonicity. 40+yr, experts, CA rigidity fully
   mapped (Wu 2019 reduced Robin to ~5,565 CA2 numbers), NO proof strategy found. Same obstruction as the analytic side.
   Caveat (honest): Cramér obstruction is monotonicity-specific; doesn't formally rule out ALL combinatorial lanes,
   but no live lane exists. Kept nugget: Caveney–Nicolas–Sondow "RH ⟺ 4 is the only composite extraordinary number"
   (cleanest reframe; equivalent-hard). Verdict: reframing to the multiplicative/elementary side re-clothes the zeros.
PATTERN NOW STRONG: two fresh fronts investigated (Λ program §LAMBDA_PROGRAM_SURVEY; elementary Front A) — BOTH reduce
to the prime-counting error / zeros. Consistent with this repo's own RH_DIFFICULTY_CONSERVATION_MAP: difficulty is
CONSERVED across reframings. Working hypothesis for the hunt: reframing is not the lever; a crack needs a genuinely new
idea about the prime-counting error/zeros themselves (or a structure outside all known framings). Front B (quasicrystal)
= last genuinely different-flavored swing (attacks zeros' GEOMETRY not arithmetic); low prior (Jensen/LP cousin ruled
out) — orient before investing.

### Front B (quasicrystal / Lee-Yang) — CLOSED as a bypass (2026-08-10, oriented)
DEAD as a route. (1) Kurasov–Sarnak (2020, stable/Lee-Yang polynomials) and the recent Fourier-quasicrystal
breakthroughs produce quasicrystals from FINITE exponential polynomials; their own statement: "among these Fourier
quasicrystals there are NONE corresponding to the zeros of the zeta function." The recent machinery stops short of ζ
(infinite Dirichlet series). (2) "Zeta zeros form a 1D quasicrystal" is CONDITIONAL on RH (a consequence, not a route);
unconditionally the zeros are only a crystalline measure via the explicit formula. (3) Dyson's actual proposal =
classify the infinite 1D quasicrystals — wide open, as hard as ever, tools don't reach ζ. (4) The one claimed proof
(Shaughnessy 2410.03673) = 12 revisions since 9/2024, unvetted = standard flawed-claimed-proof signature. Its cousin
(Jensen/Laguerre–Pólya) already ruled out here ("TP2 fails at corner"). Verdict: re-clothes RH as an equally-hard
open classification; no leverage.

## MAP STATUS (2026-08-10): ALL accessible fronts investigated & closed as bypasses.
23 analytic dead-ends + Λ program + Front A (elementary) + Front B (quasicrystal) — EVERY one reduces to the same
obstruction (prime-counting error / the zeros). RH_DIFFICULTY_CONSERVATION_MAP confirmed empirically: difficulty is
conserved across ALL known reframings. HONEST CONCLUSION: reframing is not the lever. A crack requires either (a) deep
new mathematics within a framing (F₁/condensed-Hodge frontier — beyond this setup), or (b) an unexpected EMPIRICAL
pattern discovered by computation that seeds a new idea (long shot, specialist-dominated). No reframe remains untried.
