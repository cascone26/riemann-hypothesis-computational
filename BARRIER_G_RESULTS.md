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

## G3 — matched-local-data pair: BUILT at finite truncation + sketch for the full theorem

**Built (same night, second pass):** `operators/barrier_g3_matched_pair.py` — two systems
discretized by the same deterministic greedy rule, targets identical below Xs=1000, system B's
target adds the zero-at-ρ term (ρ=0.75+2i) only above Xs. **Output:** 168 g-primes below the split,
**bit-identical** (max |p_A−p_B| = 0); above it ψ_B−x carries the prescribed x^{0.75} oscillation
(residual ≤ 8.1 after removing it) while ψ_A−x stays within O(log X) (≤ 9.7). A certificate reading
local data below Xs provably returns the same verdict on both. Also
`operators/barrier_g2b_discretize.py`: the witness system alone at X=10⁶ (77,941 g-primes,
oscillation amplitude 10,188 vs residual ≤ 12.2 < log X = 13.8, flat control system as null).
The discretization machinery is **handled on our instrument**, not just cited.

**THEOREM (G3, PROVEN — elementary; third pass 2026-08-07).** *For every X ≥ 16 there exist two
discrete Beurling g-prime systems P_A, P_B with P_A ∩ [1,X] = P_B ∩ [1,X] (identical local data)
such that ζ_A is zero-free in Re s > 1/2 while ζ_B has zeros exactly at ρ = 3/4 ± 2i (simple).
Consequently no certificate reading g-prime data up to any finite height decides the RH-analogue.*

Proof. (1) *Construction.* Targets Ψ_A(x) = x − x₀, Ψ_B = Ψ_A below X and
Ψ_B(x) = (x−x₀) − 2Re[(x^ρ − X^ρ)/ρ] above X (x₀ = 16; density Ψ'_B = 1 − 2x^{σ₀−1}cos(t₀ log x) ≥
1 − 2·16^{−1/4} = 0 for x ≥ 16, and > 0 for x > 16 off a measure-zero set). Greedy discretization:
p_{n+1} solves Ψ(p_{n+1}) = Σ_{i≤n} log p_i + log p_{n+1}; the solution exists and increases since Ψ
is continuous increasing and unbounded. Determinism + identical targets on [x₀, X] ⟹ identical
primes below X (verified bit-identical numerically at Xs=1000).
(2) *Error bound.* At each placement, ψ_P(p_n) = Ψ(p_n) exactly; between placements ψ_P is constant
while Ψ rises by ≤ log p_{n+1} ⟹ E(x) := ψ_P(x) − Ψ(x) satisfies |E(x)| ≤ log(2x) for all x.
(3) *Transfer to the zeta.* For the discrete system, −ζ'_P/ζ_P(s) = ∫ x^{−s} dψ_P^{(Λ)} where
ψ^{(Λ)} adds prime-power terms; the powers contribute Σ_p Σ_{k≥2} log p · p^{−ks}, absolutely
convergent and analytic in Re s > 1/2 (g-primes have ~linear density ⟹ Σ p^{−2σ} < ∞ for σ > 1/2).
For the prime part, Stieltjes parts: ∫ x^{−s} dψ_P = ∫ x^{−s} dΨ + [x^{−s}E(x)] + s∫ E(x)x^{−s−1}dx;
since |E| ≤ log 2x the last integral converges absolutely and is analytic in Re s > 0. And
∫ x^{−s}dΨ_B = 1/(s−1) − 1/(s−ρ) − 1/(s−ρ̄) + (entire in Re s > 0 corrections from finite lower
limit), ∫ x^{−s}dΨ_A = 1/(s−1) + (entire), by direct computation.
(4) *Conclusion.* In Re s > 1/2: −ζ'_A/ζ_A is analytic except the pole at s=1 ⟹ ζ_A has no zeros
there. −ζ'_B/ζ_B = −1/(s−ρ) − 1/(s−ρ̄) + analytic-except-s=1 ⟹ simple poles with residue −(−1)=+1
at ρ, ρ̄ ⟹ ζ_B has simple zeros exactly at ρ and ρ̄. ∎

(Sign bookkeeping: a zero of ζ of multiplicity m gives −ζ'/ζ a pole of residue −m... fix orientation:
ζ'/ζ has residue m at a zero, so −ζ'/ζ has residue −m; our computed residue of −ζ'_B/ζ_B at ρ is −1
(the −1/(s−ρ) term), consistent with a SIMPLE ZERO m=1. Checked both directions; the displayed
formula −ζ'/ζ = Σ Λ n^{−s} = 1/(s−1) − Σ_ρ 1/(s−ρ) + ... matches the classical explicit-formula sign.)

Same machinery, one line, upgrades G2: each family member can be discretized with its zeros exactly
as prescribed (zero-free good members via Ψ = x−x₀ targets shifted, the bad member via its ρ-term),
so the averaging-wall family exists as genuine discrete systems, not only zero configurations.

Status ladder note: this is our own elementary proof, self-checked including the sign audit —
own-math-unreviewed until an independent skeptic pass (standing table discipline, cf. §791); flagged
for adversarial verification next session. The numerics above (bit-identical primes, oscillation
carried, control flat) are its instrument-level confirmation.

Original sketch (kept for provenance):

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

## G4 — CLOSED for the integrality class (Hamburger assembly, 2026-08-07 second pass); general class open

**Theorem (assembled; candidate-novel as an observation — lit-gate found no published Beurling
application of Hamburger).** Let (P,N) be a Beurling system whose g-integers lie in ℕ (the
Pintz–Révész integrality condition — the same class their zero-density theorem lives in), with
counting multiplicities G(ν), such that ζ_P(s) = Σ G(ν)ν^{−s} converges absolutely for Re s > 1,
(s−1)^m ζ_P is entire of finite order, and ζ_P satisfies Riemann's exact functional equation
π^{−s/2}Γ(s/2)ζ_P(s) = π^{−(1−s)/2}Γ((1−s)/2)ζ_P(1−s). Then ζ_P = c·ζ **(Hamburger 1921, Math. Z.
10, 240–254; hypotheses verified against source 2026-08-07)**; a Beurling system has G(1)=1 (empty
product), forcing c=1, so ζ_P = ζ and, by uniqueness of Dirichlet coefficients and of the Euler
product, P = the rational primes.

**Consequence for the map:** in the integrality class, the property-set {Euler product, Riemann FE}
admits NO RH-violating witness unless RH itself is false — the empty witness column of surviving-row
1 is a THEOREM there, not an observation. The rigidity is total: the only object in the class is ζ.

**Honest scope:** (i) this is for the EXACT Riemann FE — Knopp (Invent. Math. 117 (1994) 361–372)
shows other functional equations can have infinite-dimensional solution spaces, so the barrier
statement does not transfer to arbitrary FEs by this route (though those solutions generically lack
Euler products — the EP hypothesis does separate work there, unexplored); (ii) g-integers ⊂ ℕ is a
real restriction — the general Beurling class (real g-integers) remains OPEN, and is genuinely
harder because no FE is even defined/known there (Révész: "we lack a functional equation");
(iii) assembly of two classical theorems — the value is the placement in the barrier lattice, and
the observation is easy once stated (likely folklore-adjacent; claim it as organization, not depth).

**The residual open target (the real G4):** formulate + prove an FE-rigidity statement for general
Beurling systems (real g-integers) — or construct the first Beurling system with any genuine FE.
Either result would be new mathematics; the map points all remaining weight here and at B18.

Original open statement (kept for provenance):

**Statement to attack:** let (P,N) be a Beurling system whose ζ_P (i) has EP over g-primes (automatic),
(ii) continues meromorphically to ℂ with finite order, (iii) satisfies a Riemann-type FE
ξ_P(s) = ξ_P(1−s) with ξ_P(s) = G(s)ζ_P(s), G a gamma-type factor. Prove: ζ_P has no zeros with
Re s > 1/2 + δ for some δ ≥ 0 — or even the rigidity statement that such systems are exactly the
arithmetic ones (rescaled rational primes / Dirichlet-type). Either result converts surviving-row 1
(EP+FE has an empty witness column) from an observation into a theorem in the one setting where FE
is the rare property. This is a Selberg-class-conjecture fragment; no known Beurling system has a
proven FE (fractal-membrane duality is conjectural — Hilberdink–Lapidus 2006). Deepest of the four;
stays open with this precise target on the queue.

## G3 — ADVERSARIAL SKEPTIC PASS (2026-08-07, anchor-thread session, independent of the barrier session that proved it)

VERDICT: **CONFIRMED — the theorem is correct**, with ONE fixable imprecision in the written error
bound. Checked all three flagged points + hunted for a fatal hole:

1. **Greedy error bound |E(x)| ≤ log(2x) — OVERSTATED near x₀, conclusion survives (the one real fix).**
   E(x) = Ψ(p_n) − Ψ(x) on [p_n,p_{n+1}), so |E| ≤ log p_{n+1}, and log p_{n+1} ≤ log(2x) ⟺
   p_{n+1} ≤ 2p_n (gap-doubling). For system A (density ≡1) and the TAIL of B (density→1) this holds
   (gaps ~ log x ≪ x). BUT for system B near x₀=16 the density Ψ_B' = 1 − 2x^{−1/4}cos(2 log x) DIPS
   TO 0, so Ψ_B rises slowly and a greedy gap there can have p_{n+1} ≫ 2p_n — the stated ∀-x bound
   FAILS on a compact initial region. FIX (standard, non-fatal): split E = E·1_{[16,X₁]} + E·1_{tail};
   the compact piece contributes an ENTIRE function to ∫x^{−s}dE (bounded integrand, finite interval),
   the tail obeys |E|≤log2x. Analyticity of the error integral in Re s>1/2 — the only thing the proof
   actually needs — is preserved. The written "|E|≤log(2x) for all x" should read "for x≥X₁, plus a
   bounded compact piece." Imprecision in the writeup, not a hole in the result.
2. **Prime-power analyticity in Re s>1/2 — CORRECT.** g-primes have ~linear density ⟹ Σp^{−2σ}<∞ for
   σ>1/2 ⟹ Σ_p Σ_{k≥2} log p·p^{−ks} absolutely convergent/analytic there. ✓.
3. **Residue sign audit — CORRECT, re-derived independently.** ∫_X^∞ x^{−s}dΨ_B picks up
   −2Re∫x^{ρ−1−s}dx = −X^{ρ−s}/(s−ρ) − X^{ρ̄−s}/(s−ρ̄) = −1/(s−ρ) − 1/(s−ρ̄) + entire (the X^{·−s}−1
   parts are removable at the poles). Residue −1 of −ζ'_B/ζ_B ⟺ residue +1 of ζ'_B/ζ_B ⟺ SIMPLE ZERO
   (m=1). Matches the classical −ζ'/ζ = 1/(s−1) − Σ_ρ 1/(s−ρ) + … sign. ✓.
4. **"Zeros EXACTLY at ρ,ρ̄" — CORRECT.** In Re s>1/2, −ζ'_B/ζ_B has poles only at s=1 (from 1/(s−1)),
   ρ, ρ̄; all other pieces (prime-powers, error integral) analytic there ⟹ no other zeros. ✓.
5. **Both members are legitimate Beurling systems.** Ψ_B'≥0 for x≥16 (⟹ greedy primes exist,
   increasing); identical targets on [x₀,X] ⟹ identical primes below X (determinism). ✓.
6. **Corollary sound:** a matched pair with opposite verdicts for every X ⟹ no finite-local-prime-data
   certificate decides the RH-analogue. ✓.

NET: G3 stands as a proven theorem. Required edit: restate the error bound with the compact/tail split
(one paragraph). Skeptic found no witness-, premise-, circularity-, or sign-error killing the result —
the independent pass the §791/§805 discipline demanded. (Reviewer: anchor-thread session, not the
barrier session that authored it.)

## §812 — the missing lattice witness (EP + ℕ, no FE): the row is STRUCTURALLY PROTECTED, and here is the proof

Target was: build an RH-violating witness with Euler product + g-integers ⊂ ℕ + Axiom A counting but
no functional equation (the lattice's one "OPEN-leaning-yes" row). The build attempt produced the
opposite — a rigidity proposition showing why every reachable construction FAILS, which flips the
row's meaning: integrality alone (without FE) already pins strip-zeros to ζ's own.

**Proposition (surgery rigidity; own-math, elementary — flagged for skeptic pass).** Let the
g-prime multiset P be obtained from the rational primes by deleting a set D and adding a multiset A
of integers ≥ 2, with Σ_{n∈D∪A} n^{−σ₀} < ∞ for some σ₀ < 1 (sparse surgery; σ₀ < 1/2 for
strip-wide conclusions). Then
  ζ_P(s) = ζ(s) · Π_{p∈D}(1−p^{−s}) · Π_{n∈A}(1−n^{−s})^{−1},
and the two products converge absolutely and are analytic in Re s > σ₀ with zeros and poles ONLY on
the line Re s = 0 (each factor 1−n^{−s} vanishes exactly at s = 2πik/log n). Hence **in the strip
Re s > max(σ₀, 0), ζ_P and ζ have IDENTICAL zero sets.** In particular no sparse surgery on ℕ can
create an off-line zero — the system violates the RH-analogue iff ζ violates RH.

Proof. The factorization is immediate from the Euler product (absolute convergence rearrangement in
Re s > 1, then analytic continuation of the finite-or-absolutely-convergent products to Re s > σ₀).
An absolutely convergent product Π(1−a_n(s)) with Σ|a_n(s)| < ∞ vanishes only where a factor
vanishes; 1−n^{−s} = 0 ⟺ n^{−s} = 1 ⟺ s = 2πik/log n, Re s = 0. ∎

**Why dense surgery can't rescue a witness:** Axiom A with positive density a > 0 forces
Σ_{p∈D} 1/p < ∞ (else the g-integer density a = lim N(x)/x picks up the factor Π_D(1−1/p) → 0),
i.e. deletions must be summable-sparse; additions with Σ 1/n = ∞ inflate N(x) beyond ax + O(x^θ).
So the Axiom-A-compatible surgeries are exactly the sparse ones the proposition covers (the gap
between "Σ1/n < ∞" and "Σ n^{−σ₀} < ∞ for σ₀ < 1/2" is real and stated honestly — surgeries in
that gap have products analytic only in a smaller half-plane, protecting a correspondingly smaller
strip region).

**Lattice consequence (row rewritten):** EP + ℕ + A (no FE) has NO reachable witness: every
multiplicative surgery on ℕ inherits ζ's strip-zeros verbatim, so producing a witness in this class
is AT LEAST as hard as exhibiting an off-line zero of ζ itself. Combined with Hamburger (EP+FE+ℕ ⟹
exactly ζ): **integrality is rigid with or without the FE — the natural numbers themselves carry
the RH-hardness.** This also explains, structurally, why every literature witness (W3–W7) uses real
g-primes: leaving ℕ is not a convenience, it is the only door.

Honest scope: "reachable" = surgery constructions (delete/add/multiplicity on integer support).
A non-surgery integer-supported system is one whose g-prime multiset differs infinitely from the
primes with only Σ n^{−1} < ∞ control — the proposition's σ₀-gap; no construction principle for a
witness there is known, and none is claimed impossible. REPORT.

**§814 Lean upgrade:** the proposition's zero-location step is now KERNEL-CHECKED —
`RiemannLean/RiemannLean/Cand_surgery_factor_zero.lean`: `surgery_factor_zero_re_eq_zero`
(1 − n^(−s) = 0 ⟹ Re s = 0, real n ≥ 2) and its contrapositive, both passing `#print axioms` on
{propext, Classical.choice, Quot.sound} only, first compile. The factorization/absolute-convergence
step of the proposition remains own-math (skeptic-pass queue with G3).
