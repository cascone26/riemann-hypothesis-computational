# Barrier gaps G1–G4 — execution results (2026-08-07, barrier session, LEDGER §803)

Status ladder per row: **built** (ran here, output cited) / **assembled** (pure citation-assembly,
sources named) / **sketch** (proof outline, load-bearing step verified) / **open** (stated precisely).

## G2 — averaging-wall witness family: BUILT (continuous level)

`operators/barrier_g2_averaging_witness.py`. Family of K=12 continuous Beurling-type zero
configurations: 11 with all zeros on the critical line (low-lying ordinates from ~1), 1 with an
off-line pair ρ = 0.95 ± 2.0i (DMV territory: σ near 1, low ordinate) and its own on-line zeros
starting at 10.9. Gaussian test h(r)=exp(−r²/w²), Weil zero-side form W_j = Σ h(γ).

**Output (verified, w = 0.42 through ≥1.06):** W_bad < 0 (individual violates positivity),
mean_j W_j > 0 (family average positivity holds), outlier at **0.61σ** — Chebyshev bound 2.69 vs
family rate 0.083, so **no variance-strength (large-sieve-type) bound can exclude the single
violator.** §62's logic-ceiling ("forcing every χ would prove GRH") is now realized as an explicit
family: mean + variance information provably never reaches the individual.

Geometry lesson (deposit): the off-line pair's term dominates its member's form only when its
ordinate t₀ sits BELOW the member's first on-line ordinate (term ratio exp((γ₁²−t₀²+δ²)/w²));
first attempt with t₀=30 above γ₁≈14 found no admissible width — the fix was the ordering, not
the width scan.

Tier: REPORT (continuous zero configurations; discretization to genuine g-prime systems = BDR
arXiv:2309.01567 Thm 1.2, cited not executed).

## G3 — matched-local-data pair: SKETCH (load-bearing step verified)

**Claim:** for every X there exist two Beurling systems with identical g-primes on [1,X] and
opposite RH-analogue verdicts. Hence no certificate reading local prime data to ANY finite height
decides the RH-analogue — upgrading B4/B8 from premise-false kills to witness-barred.

Sketch: take F_A(x) = li(x) (target of an RH-true continuous system) and
F_B(x) = li(x) − [li(x^ρ) + li(x^ρ̄)]·1_{x>X} (smoothed at X), ρ = σ₀ ± it₀. Discretize both by BDR
Thm 1.2 with the SAME random draw on [1,X] (the selection algorithm processes [1,X] identically for
identical targets) ⇒ g-primes literally agree below X.

**Load-bearing step (verified here):** the zero at ρ survives restriction of the li(x^ρ) term to
x > X, because the log-singularity of ∫ x^{−s} d li(x^ρ) at s=ρ comes from the divergence at ∞,
not the lower limit: ∫_X^∞ t^{ρ−1−s} dt/log t = −log(s−ρ) + (entire near ρ) for any X. The
omitted [1,X] chunk is entire ⇒ ζ_B still has its zero at ρ; ζ_A is zero-free off-line.
Positivity of the prime density needs 2x^{σ₀−1}cos(t₀ log x) ≤ 1, i.e. X ≥ 2^{1/(1−σ₀)}
(X ≥ 16 at σ₀ = 0.75) — or BDR's M·li(x^δ) compensator in general. Tier: REPORT (sketch).

## G1 — zero-density-strength + off-line zero: RESOLVED — ALREADY IN THE LITERATURE (stronger form)

The lit-gate paid off exactly as designed. The citation-assembly works (BDR Thm 3.2 prescribes a
zero at any σ₀ ∈ (1/2,1) in a system satisfying Axiom A with θ=β; Révész's Carlson-type theorem
arXiv:2209.01689 Thm 1 applies under Axiom A ALONE, no FE — N(σ,T) ≤ 1000(A+κ)⁴/((1−θ)³(1−σ)⁴) ·
T^{12(1−σ)/(1−θ)} log⁵T; a finite prescribed zero set is trivially consistent). **But the
conclusion is published, in stronger form:** Révész 2209.01689 §1.2 explicitly cites
**Broucke–Debruyne, Acta Arith. 207 (2023) (arXiv:2211.08716)**: a Beurling system under Axiom A
with **Ω(T^{c(1−σ)}) zeros** in [σ,1]×[−iT,iT] for σ arbitrarily close to 1 — with the remark that
"Carlson-type density estimates might be seen in general as the best possible." Not one off-line
zero: power-of-T many, saturating the density bound.

**Disposition:** B5's witness column is FILLED from the literature (new W7); the barrier
"zero-density-strength information cannot force RH" is DEAD-proven, not ours to claim. Exactly the
rh-roundtable-literature-gate failure mode avoided. Bonus banked from the same pass: Révész notes
the zero-density theorems succeed WITHOUT any functional equation (the Selberg-class treatment
needs one) — reinforcing that counting/density machinery is FE-blind, which is why it walks into
witness territory.

## G4 — FE + EP ⟹ RH-analogue (narrow Beurling class): OPEN, stated precisely

**Statement to attack:** let (P,N) be a Beurling system whose ζ_P (i) has EP over g-primes (automatic),
(ii) continues meromorphically to ℂ with finite order, (iii) satisfies a Riemann-type FE
ξ_P(s) = ξ_P(1−s) with ξ_P(s) = G(s)ζ_P(s), G a gamma-type factor. Prove: ζ_P has no zeros with
Re s > 1/2 + δ for some δ ≥ 0 — or even the rigidity statement that such systems are exactly the
arithmetic ones (rescaled rational primes / Dirichlet-type). Either result converts surviving-row 1
(EP+FE has an empty witness column) from an observation into a theorem in the one setting where FE
is the rare property. This is a Selberg-class-conjecture fragment; no known Beurling system has a
proven FE (fractal-membrane duality is conjectural — Hilberdink–Lapidus 2006). Deepest of the four;
stays open with this precise target on the queue.
