# RH ATTACK QUEUE — the anti-stall mechanism

**Purpose (Jacob, 2026-07-18):** never hit a "what do I do next?" pause. This is a LIVE,
ranked queue of genuinely-NEW attack methods — each attacks the actual wall (Weil positivity
A+B≥P) or opens an orthogonal front, NOT another statistical confirmation. Self-drive: pop the
top unblocked item, RUN it, record result in LEDGER, mark it here, refill the queue. Always
keep ≥5 live items. Read this at session start with LEDGER.md.

**Prime lesson (07-18):** statistics on the zeros can't crack RH (every probe confirms; the
forcing is analytic/global, not statistical). So NEW methods must attack POSITIVITY STRUCTURE
directly, or find a genuinely orthogonal handle. Rank by: (P) proof-relevance, (N) novelty
(not tried in this project), (R) runnable-now.

---

## ACTIVE QUEUE (ranked)

### A1. Certified Cholesky/LDL^T STRUCTURE mining  [P:high N:high R:now]  <-- RUNNING
The Groskin tool certifies τ⪰0 via interval LDL^T. Nobody has looked at the PIVOT VALUES /
Cholesky factor as objects. Extract d_i(N,c) and L_ij at high precision; hunt for a closed form
or manifest-positivity pattern (each pivot = sum of squares? a recognizable recurrence?). A
clean structural law for why the pivots stay positive uniformly in N = a CONSTRUCTIVE positivity
proof seed. This is the most direct attack on the wall. Script: `operators/ldlt_structure_mine.py`.

### A2. Identify the critical (near-kernel) eigenvector in CLOSED FORM  [P:high N:med R:now]
weil_cancellation_analysis showed the min-eigenvector converges to a fixed limit function g_lim
as N→∞. Use the high-precision cutoff-free Groskin τ to nail g_lim's coefficients to 100+ digits,
then symbolic-recognize them (OEIS / inverse-symbolic / ansatz fits: constant? prolate? a known
special function?). If g_lim is identified, prove B(g_lim) ≥ |A(g_lim)| + P(g_lim) analytically →
closes the wall on the tightest direction. Script: `operators/critical_eigvec_id.py`.

### A3. Symbolic regression on the A/B/P decomposition  [P:med N:high R:now]
Feed the analytic A(e_n), B[n,m], P[n,m] matrix entries (already have digamma_poisson_decomp) to
symbolic regression (PySR if installable, else ansatz-fit) to find the exact structural inequality
that makes A+B−P PSD. New tool on an old object; may surface a relation the hand-analysis missed.

### A4. de Bruijn-Newman BACKWARD heat-flow margin  [P:high N:high R:now]
RH ⟺ Λ=0, Λ≥0 proven. Under backward flow (t<0) ξ_t zeros LEAVE the line if RH is fragile. Compute
backward-flow zero trajectories from real zeros, measure how far back (t→Λ⁻) they stay real / the
collision margin. A DYNAMICAL probe (not static stats): quantifies how much "room" RH has, and the
margin's N-scaling could show whether Λ is pinned exactly at 0. Extend `operators/debruijn_newman.py`.

### A5. Lee-Yang positivity transfer, made concrete  [P:high N:high R:med]
Physics dive: ξ has Lee-Yang's reflection symmetry (funct. eq.) but lacks the coefficient POSITIVITY
Lee-Yang uses. Build the Lee-Yang-style positivity for a FINITE truncation of ξ (as a partition
function w/ product structure) and localize EXACTLY where positivity breaks — that break IS the RH
content, now pinned to a concrete finite algebraic obstruction. New framing of the wall.

### A6. GRH cross-L-function universality of certified-PD  [P:med N:high R:now]
Does the certified τ⪰0 hold for Dirichlet L-functions (GRH) with the analogous construction? If the
positivity mechanism is universal across the family, that's a strong structural constraint on any
proof (it can't use anything ζ-specific). Adapt the Groskin build to a small conductor L-function.

### A7. Prime-adapted / wavelet basis where P is manifestly dominated  [P:med N:med R:now]
All prior work used sinc/PSWF bases. Try a basis adapted to the prime points {log p} (a Gabor/wavelet
frame) where the prime form P becomes near-diagonal or manifestly ≤ A+B. If such a basis diagonalizes
the obstruction, positivity may become visible. New basis = new chance.

### A8. Extract a rigorous N-UNIFORM lower bound on the near-kernel pivots  [P:high N:med R:now]
Not the smallest eigenvalue (→0), but: is the k-th smallest CERTIFIED pivot bounded below by an
explicit positive f(N) whose form is N-independent-enough to extend to ∞? Interval arithmetic gives
rigorous bounds the old float64 work couldn't. If a uniform structural floor exists on the RIGHT
subspace, that's the proof. (Pairs with A1.)

---

## RUN THIS SESSION (2026-07-18) — outcomes

- **A1 (LDL^T pivot mining):** all pivots certified >0 (re-confirms PD), BUT natural-order pivots shrink
  erratically with N and the near-kernel spreads across the ordering → raw pivots are NOT a clean invariant.
  `results/ldlt_structure_mine.json`. Outcome: no structural handle from natural-order pivots.
- **A2 (identify critical eigenvector):** ILL-POSED — near-kernel is 23-dim DEGENERATE, so the min-eigval
  eigenvector is an arbitrary vector in that subspace (coeffs don't converge across N). No single critical
  direction exists. `results/critical_eigvec_id.json`. Spawned A9.
- **A9 (low-freq Schur decomposition):** τ⪰0 ⟺ Hi⪰0 (high-freq) AND S⪰0 (low-freq Schur complement). Real
  findings: (i) S CONVERGES to a fixed finite S∞ as N→∞ (good — localizes the wall to a finite object);
  (ii) Hi is PD with a fixed N-independent gap (~3e-8 at M=11 → 5e-5 at M=15) but diag-dominance margin is
  NEGATIVE (positivity NOT manifest — Gershgorin fails); (iii) S∞ has a ~19-20 dim KERNEL (out of 23-31),
  so the degeneracy SURVIVES the Schur complement. `results/low_freq_schur.json`. Outcome: clean
  reformulation, wall relocated to fixed finite S∞⪰0 with ~19-dim kernel, but NOT dissolved.

**META-LESSON (2026-07-18, important):** A1/A2/A9 are three linear-algebra views (pivots / eigenbasis /
Schur), and ALL THREE hit the SAME 23-dim positivity-boundary content in different coordinates — exactly
what OBSTRUCTION_MAP predicts. **A fourth linear-algebra reformulation will hit it again. Next attacks must
be a DIFFERENT KIND of method** (dynamical/heat-flow A4, physical-positivity-transfer A5, cross-family A6,
or a non-orthonormal prime-adapted frame A7), not another matrix reshuffle. Prioritize A4 (orthogonal:
dynamical) next.

## NEW LEADS — run this session, OUTCOMES
### A10. High-freq block Hi as H-matrix?  DONE → NO. (`operators/hi_block_hmatrix.py`)
Hi is NOT an H-matrix (min(M(Hi)⁻¹·1)<0, only ~half rows diag-dominant) at every N,M. So Hi⪰0 is NOT
provable by diagonal scaling; the high-freq block's positivity is genuinely nonlocal too. Wall does NOT
reduce to S∞ alone. Negative.
### A11. Kernel dim = arithmetic count?  DONE → NO clean identity. (`operators/kernel_dim_arith.py`)
Near-kernel dim grows ~linearly in c (slope ≈1.91, ~2c − slow correction), NOT as π(c) / #prime-powers
(both ~c/log c, far smaller). So the kernel fills a whole frequency BAND; primes are sparse within it →
reinforces "positivity is global, not one-dimension-per-prime." No arithmetic formula.

### A12. SUSY / superpotential factorization τ = A†A  [P:high N:high R:med] <-- NEW TOP ORTHOGONAL FRONT
From Jacob's "antimatter/reverses/counter" intuition (2026-07-18), grounded in real lit. In SUSY-QM,
H = A†A is MANIFESTLY ≥0 (positivity built in, not proven) — exactly what the wall needs and what 5 linear-
algebra reshuffles (A1/A2/A9/A10/A11) provably CAN'T give. We tried numerical Cholesky (no structure); we
did NOT try a STRUCTURED first-order factorization A = (shift/difference operator) + W (superpotential from
the primes). Goal: find A with arithmetic meaning s.t. τ=A†A → constructive positivity proof.
- Anchors: arXiv:2301.05360 (2023, High-Order SUSY-QM + XP model: SUSY partner H's whose ZERO-energy mode
  sits on the critical line; superpotential W=-d/dζ ln ψ₀, |x|^{-s} seed, factorization energy 0). Spector
  1990 (CMP): Möbius μ(n) = Witten index (-1)^F = SUSY signed count; 1/ζ = fermionic partition function.
- Coherence: Lee-Yang (§6) = the REQUIREMENT (reflection+positivity); SUSY = the MECHANISM (A†A factorization);
  functional eq. ξ(s)=ξ(1-s) = the reflection/SUSY symmetry. Fits the whole picture.
- First runnable probe: does the (regularized, near-kernel-projected-out) Cholesky factor of τ approximate a
  FIRST-ORDER operator (banded: diagonal + one off-diagonal = discretized d/dx + W)? If L is near-first-order,
  read off the superpotential W_n and check if it's arithmetic (built from log p / prime positions). If L is
  dense/structureless, SUSY factorization in THIS basis is out (try prime-adapted basis A7 first).
- HONEST: moonshot (a clean factorization ≈ the proof). But it's the RIGHT kind, aimed at the actual wall.
- Möbius/Witten sub-lead: signed-count RH criteria (Mertens M(x), Pólya) are numerically weak (both
  conjectures DISPROVEN at astronomical x) — structural insight, not a computable lever. Don't over-invest.

### A12 SUSY factorization — DONE (probe). (`operators/susy_factorization_probe.py`, `/tmp/susy_W.py`)
BEST structural find of the session: τ^½ (a valid SUSY A=A†) is QUASI-LOCAL — ~98% band-energy at |d|≤1,
off-diag decay 0.13/0.05/0.04, Toeplitz score ~1.3 (n-varying diagonal = a superpotential W(n)). It cleanly
separates near-kernel (W_n~0.02 on low modes) from bulk (W_n~1.5). So the Weil operator HAS an approximately-
local square root — a genuine SUSY-like signature, NOT the pure nonlocality the other 5 attacks showed.
BUT: W_n is jagged/not-cleanly-arithmetic in the frequency basis (W_6=1.62,W_8=0.13,W_10=0.69…), and the
nonlocal tail (mode-0 couples out to distance 22 = the prime/RH content) is spread, not sharp. No clean
factorization → no crack. Real structure, honest partial-positive.
### A13 (spawned from A12) — PARTLY DONE. (`operators/susy_basis_A13.py`)
Tested arch (A+B) eigenbasis: superpotential jaggedness 0.94→0.58 (~40% smoother) AND a PARITY structure
appears (τ^½ couples even-separated modes ~0.13/0.065/0.044, odd ~0.001 — even/odd sectors ~decouple). So W
IS basis-dependent and cleans in a natural basis, but 0.58 still not clean (want ≪0.1). NEXT (careful build,
next session): position/log-line basis where primes localize at x=±log p — build τ in position space (arch =
smooth convolution kernel Ω̌(x)=−e^{−|x|/2}/(1−e^{−2|x|}), prime = point masses at log p, boundary = rank-2),
take √, read W(x). If W(x) is clean/arithmetic there → explicit SUSY factorization. Highest-life thread.

## SESSION CONCLUSION (2026-07-18): LINEAR-ALGEBRA FAMILY EXHAUSTED
Five attacks — A1 (pivots), A2 (eigenvector), A9 (Schur), A10 (H-matrix), A11 (kernel-count) — ALL confirm
the τ positivity is genuinely NONLOCAL: it does not become manifest under any pivot ordering, eigenbasis,
block/Schur decomposition, diagonal scaling, or kernel-counting. This is strong, multi-angle re-confirmation
of OBSTRUCTION_MAP's thesis from new directions. **STOP doing linear-algebra reshuffles of τ.** The only
attacks that could escape this are a DIFFERENT KIND:
  - **A4 (dynamical):** de Bruijn-Newman heat-flow zero dynamics (orthogonal; but mostly re-derives Λ≈0 — low).
  - **A6 (cross-L universality) — DONE 2026-07-18:** τ_χ for L(s,χ₄) sits EXACTLY on PSD boundary (min eig
    ~−3e-14) like ζ, despite no pole; knife-edge at the correct conductor const log(4/π). POSITIVITY IS
    UNIVERSAL across L-functions → no proof can use ζ-specific structure. `operators/cross_L_A6.py`.
  - **A5 (Lee-Yang positivity transfer):** localize where ξ's missing coefficient-positivity breaks, concretely.
  - **A7 (prime-adapted frame):** still linear-algebra-ish; lower priority given the exhaustion above.

## DONE / MOVED TO LEDGER
- Certified τ⪰0 on c×N grid (Groskin) — LEDGER §4/§4b.
- Jensen quartic discriminant bridge 27Δ=4I³−J² — LEDGER Track 3.
- Structure factor, HKO moments (RS evaluator), Lehmer gaps, FHK max — LEDGER §6/§7 (all consistency, no crack).
- Physics deep-dive (Lee-Yang=wall) — `research/RH_PHYSICS_DEEPDIVE_2026-07-18.md`.

## REFILL RULE
When queue < 5 active items, generate more by asking: (1) what NEW representation makes positivity
manifest? (2) what rigorous tool (interval arith, exact symbolic, CAS) can now do what float64 couldn't?
(3) what orthogonal object (L-function family, heat flow, dual cert) constrains the proof? Never sit idle.

## NEW METHODS (2026-07-18, from solved-analog mining) — the better generative engine
**PROCESS (use this to generate methods, not mine our own failures):** for each SOLVED problem
isomorphic to RH — function-field RH (Weil/Deligne), Selberg (hyperbolic), Lee-Yang, Jensen d≤8 —
extract the EXACT mechanism that made it work, find the one ingredient missing in the number-field
case, build just that. Obstruction map = closed doors here; this = open mechanisms elsewhere + what's missing.

### A14. Toeplitz total-positivity (PF sequence) — DONE, valid, extension fails. (`operators/toeplitz_tp_test.py`)
RH ⟺ Ξ coeff sequence is Polya-Frequency ⟺ its Toeplitz matrix is Totally Positive (ASW-Edrei). DIFFERENT
structure than dead Hankel/moment. TESTED: Ξ coeff Toeplitz IS TP to n=15 (Neville all-positive); TP2 ratios
2.15→1.02 decreasing to 1 from above (holds ∀n, asymptotically tight). Payoff idea (Euler product → TP-closure)
FAILS: ξ=arch×ζ but ζ(½)<0 so factors aren't locally positive; positivity is product-emergent (arch↔prime
cancellation). KEY: 4 frameworks now (moment/operator/SUSY/TP) all show positivity is IRREDUCIBLY product-level,
never factorable into local/per-prime pieces → any local/factorization method fails for the same reason. Rules
out a whole class. Other TP handles unexhausted: variation-diminishing, sufficient conditions, the →1 asymptotic.

### A15. Function-field Hodge-index positivity transfer [P:high N:high R:hard] <-- GENUINELY DIFFERENT KIND
The one solved analog whose positivity is GEOMETRIC not analytic: Weil/Deligne prove RH-for-curves via the
Hodge index theorem — the intersection pairing on C×C is positive-definite on a subspace, forcing |α|=√q. This
is cohomological, escapes the "emergent-nonlocality/local-methods-fail" trap (A14) by construction. Transfer:
is there a computable analog intersection form for ζ whose signature (1,−,−,…) we can check? Connes-Consani's
adele-class-space is the 25-yr attempt; the specific "compute the Hodge-index signature of the analog pairing"
is the missing concrete step. Hard/deep, not a quick script — but the right KIND. Next major front.

---

## ACTIVE QUEUE — 2026-07-23 (first-prime-window / tail-lemma arc, LEDGER §184-201) — MULTI-SESSION
Supersedes the operator-positivity-family queue above for CURRENT priority (that basin's exhaustion +
the strategic pivot are logged in LEDGER/NEXT_SESSION_RH.md; kept above for provenance, not live). These
3 items are the honest forward set named in LEDGER §200's step-back audit + §201, now split across 3
concurrent tracks (see SESSIONS.md for the coordination protocol). Claim tag format:
`<-- CLAIMED by:<track> since:<ISO-ts>`; release/mark done when exhausted, per SESSIONS.md protocol.

### Q1. Tail lemma finalization (K0 certification)  [P:critical N:high R:now]  <-- CLAIMED by:tail since:2026-07-23T20:00:00Z
LEDGER §201: the k->infinity tail lemma structure is fully derived (diagonal log-growth, 1/l off-diagonal
coupling C_k=2k[g(0)+g(2a)], Gershgorin dominance for k>=11). Remaining: (1) lock C_k numerically
(coupling_asymptotic.py), (2) turn Gershgorin dominance into an explicit delta(K0) lower bound, (3)
reconcile the workflow's 6e-8 near-null cancellation vs the 1/K0 worst-case residual to pick the smallest
safe K0, (4) re-run arb_certify at that K0. This is LEDGER's own named critical path (§200 reframe: one
proven tail lemma upgrades the ALREADY-CERTIFIED Groskin c<=100 grid to true infinite-dim Weil positivity,
support up to a~2.30). Files: rigor/tail_structure.py, tail_decay.py, coupling_asymptotic.py,
arb_certify.py (extend additively, log in SESSIONS.md shared-file table first).

### Q2. Continuum-in-a feasibility (a-interval rigorous cert)  [P:high N:high R:now]  <-- CLAIMED by:continuum since:2026-07-23T20:00:00Z
LEDGER §199.6: naive a-interval Arb FAILS (entry-level over-widening, pivot blew up to +-312 over a
width-0.01 panel). Margin slope is small (0.0027) so a point-grid IS feasible (~20-40 points) if a rigorous
bound on dlambda_min/da exists. Two known-technique routes, either closes it: (A) rigorous eigenvector
enclosure (Rump/residual methods) at each grid point, or (B) first-order Taylor-model interval arithmetic
in a. Files: arb_panel.py (shows the failure mode), margin_deriv.py, arb_certify.py (point cert reference,
read-only unless coordinating with Q1/Q3 per SESSIONS.md).

### Q3. Odd sector certification + screw/CvS equivalence check  [P:med N:med R:now]  <-- CLAIMED by:oddsector since:2026-07-23T20:00:00Z
LEDGER §200 Gap 2 (flagged easiest — odd sector's single negative arch mode has 4-40x more margin than
even) + Gap 3 (honesty check: is the Suzuki-screw form's a=0.549 point-cert independent of the
Groskin/CvS grid's c=3 point, or the same object in a different basis?). Two sub-tasks: (a) extend
certification to odd modes (arb_certify.py currently does even modes [5,7,3,1] only — additive extension,
log in SESSIONS.md first), (b) cheap numeric comparison of lambda_min at matched (a,c) between the two
forms to settle Gap 3. Files: arb_certify.py (extend additively), new comparison script in
operators/track_oddsector/.

**Refill**: when a track exhausts its item, check SESSIONS.md + LEDGER's latest §-entries for the next
unclaimed relevant item before generating a new one from scratch (search-before-build applies to queue
refill too).
