# RH Research Queue — ranked, saved 2026-07-20 (riemann4, Opus)

Eight directions generated after the 07-20 meta-theorem (positivity is irreducibly global/product-emergent;
see `rh-positivity-global-metatheorem` memory, LEDGER §37, NEXT_SESSION_RH.md). Ranked in execution order.
All are non-circular and meta-theorem-compatible. "Bounded" ETAs are for the build+run; findings are open.

## Execution order (my ranking)

### #0 [NEW TOP PRIORITY 2026-07-20 → LEDGER §41] Adopt Suzuki's screw-function / Kreĭn framework
Search-before-build on #3 found the LIVE FRONTIER on our exact object: Masatoshi Suzuki, "Weil's quadratic
form via the screw function," arXiv:2606.09096 (June 2026), + his 2023 [14] (Hilbert–Pólya op in §6). The screw
function g(t) IS our A+B−P as one continuous function; RH ⟺ g positive-definite screw kernel (Kreĭn–Langer).
Gives (a) a rigorous COERCIVITY result on the localized form: λ_a=log(1/a)+μ₁−log2π+ψ(2)−1+O(a)>0 for small a
(handle (i)); (b) Thm 1.5/Cor 1.6: Weil form ≅ momentum operator 𝒟_a=i d/dx on a PRIME-WEIGHTED space ℋ(T_a),
extensions by boundary phase θ, char. fn W(a,θ;z) entire with ALL REAL zeros (Hermite–Biehler = handle (ii)),
RH ⟺ tune θ(a) so a→∞ gives z²ξ/ξ'. This RIGOROUSLY IS today's §40 conclusion (independent convergence).
NEXT: VERIFY 2606.09096 + [14] against real PDFs (current status = REPORT-tier, read via fetch); reproduce λ_a
numerically on [−a,a] (direct Thm 1.4 check + live coercivity handle); then attack Cor 1.6's θ(a)-tuning — the
sharp current target unifying handles (i)+(ii) into ONE object. Also read Connes arXiv:2602.04022 (2026).
Toolkits NEW to this project: Kreĭn screw functions, Suzuki–de Branges spaces.


### #1 [DONE 2026-07-20 → LEDGER §38] Prime-deformation cartography — treat the primes as a VARIABLE
Every prior session held the primes fixed and varied the tool. Invert it. The §37 pipeline computes A+B−P
for ANY weighted comb (validated to 30 digits vs the zero side). Deform the prime system continuously —
jitter log p by ε, thin them, swap in random Cramér primes, swap in a Diamond–Montgomery–Vorhauer "bad"
g-prime system where the RH-analog is KNOWN to fail — and measure the positivity margin as a functional over
prime-space. Question: is the true prime system a CRITICAL POINT of the margin functional (LP-duality/knife-
edge predicts yes)? Does the margin die smoothly or discontinuously as we leave the integers? WHICH feature
of the actual product buys positivity? Distinct from `beurling_adversary.py` (discrete hand-picked swaps at
finite τ_N, selectivity check) — this is a continuous path + margin functional on the certified pipeline.
DONE: ℤ is the UNIQUE conical local MAX of λ_min (peak=knife-edge 0; off-prime combs are saddles; N=6..18). Tool: `operators/prime_deformation_cartography.py`, results/prime_deformation_cartography.json. **DECISIVE FOLLOW-UP (still open):** redo with DENSITY-CONSISTENT A (co-move archimedean term with the deformed counting function) — tells real variational characterization from freeze-A artifact.

### #2 [DONE 2026-07-20 → LEDGER §39] Graph-zeta sibling — comparative anatomy vs a SOLVED world
RH for Ihara zeta of Ramanujan graphs is PROVEN (spectral gap; LPS/Margulis). Finite graphs, exact linear
algebra, no quadrature. Build the τ = A+B−P analog there (geodesic cycles = "primes", spectral/Weyl term =
archimedean), diff the mechanism against ζ term by term: in the solved world, WHAT carries the positivity ζ
lacks? DONE: mechanism diff = ζ lacks a BOUNDED SELF-ADJOINT operator (graph has the adjacency matrix, self-adjointness free; ζ's missing Hilbert–Pólya operator IS the problem). Loops to #6. Tool: `operators/graph_zeta_tau.py`, results/graph_zeta_tau.json.

### #3 [DOWNGRADED 2026-07-20 — likely folklore/frontier, see §41] χ-family average
Averaging over χ mod q is global orthogonality → escapes the meta-theorem. Prime side thins ×φ(q), arch side
grows ~log q → averaged inequality holds with margin; content = large-sieve variance on the χ-axis (parallel
to Theorem 2 §28's n-axis second moment). Candidate: positive Weil form for almost-all χ mod q, unconditional,
bandwidth ≲ log q. MANDATORY lit-check first (Iwaniec–Kowalski Ch5/7, Murty–Perelli, Sarnak–Shin–Templier) —
DOWNGRADED: Lagarias math/0404394 already links Li coeffs↔Weil functional for L-families w/ positivity+asymptotics; χ-variance is mature hard moment theory. Theorem 3 = folklore-or-frontier, not a clean win. Don't build blind. Full mechanism still in NEXT_SESSION_RH.md if revisited.

### #4 Logic lens — RH is Π₁, so independence ⟹ truth
RH is Π₁ (Lagarias elementary form: σ(n) ≤ H_n + exp(H_n)log H_n). A false Π₁ statement has a finite
counterexample = a disproof. So if RH is unprovable in PA, it is TRUE. The logic lens splits the problem:
prove it, OR prove it's independent — and the second counts as the first. Proof-mining (Kohlenbach) of
existing partial results is the workmanlike version. Low probability, nonzero, genuinely different toolset.
Reading + scoping, not a script yet.

### #5 p-adic world — Hilbert–Pólya is already TRUE at every finite place
Iwasawa Main Conjecture (proven, Mazur–Wiles): zeros of p-adic L-functions = eigenvalues of Galois actions
on Iwasawa modules — a PROVEN spectral interpretation of zeros at each place p. The missing place is ∞ (the
archimedean one) — exactly where our sign-indefinite Ω lives. Framing: RH isn't missing a spectral theory,
it's missing it at ONE place out of infinitely many proven ones; the proven cases constrain the shape of the
missing one. Reading direction.

### #6 [DONE 2026-07-20 → LEDGER §40] Physics, precise — deficiency indices computed
Berry–Keating H=xp's unfinished piece is self-adjointness of the regularized operator (Bender–Brody–Müller
2017). Self-adjointness reduces to deficiency indices + boundary conditions — computable functional analysis.
xp on the half-line has a circle of self-adjoint extensions giving the SMOOTH counting; the primes must enter
as the BOUNDARY CONDITION at the edge of the dilation group. Deficiency-index numerics on the regularized
operators — DONE: (0,∞) e.s.a.(0,0) continuous spectrum; truncated [a,b] is (1,1)=momentum-on-circle, RIGID spectrum (gap 2π/log(b/a)), can't match ζ's irregular zeros with a LOCAL boundary phase ⟹ Hilbert–Pólya operator must be NONLOCAL (periodic-orbit BC). Tool: berry_keating_deficiency.py.

### #7 Variational metaphysics, cashed out
Universality + knife-edge + LP-duality ⇒ RH behaves like the OPTIMUM of a variational principle (margin
exactly 0, every family, complementary slackness satisfied by the zeros). "What is extremized, over what space
of arithmetics?" — the computable version of this IS #1 (space = g-prime systems, functional = margin, test =
is ℤ critical). #1 is this idea's experiment; keep this as the interpretive frame for #1's output.

### #8 Far shore — reading only, flagged not romanticized
Clausen–Scholze condensed/analytic geometry as the live candidate to finally build "Spec ℤ as a curve" (Vector
D's missing object). Sarnak Möbius-disjointness / derandomization (RH = primes are computationally
pseudorandom). Long-horizon reading; no near-term script.

---
Standing constraints for ALL items (from the meta-theorem): no local/per-prime/factorization/single-positive-
test-function certificate — those are dead, cite the meta-theorem. Verify-before-ship every claimed result.
Search-before-build incl. ATTACK_QUEUE A-item OUTCOMES. Append PROCESS.md decision chain.
