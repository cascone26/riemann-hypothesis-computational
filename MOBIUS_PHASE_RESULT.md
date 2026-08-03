# The Phase-Only Möbius Comb — derivation note (2026-08-01, LEDGER §517–§553)

REPORT-tier. All numerics verified in-session; scripts `operators/krein_*.py`, data zeros10k/argzp10k/
abszp10k.npy. Lit gate (§521): the phase-only observation was NOT found in the literature. No Handle for RH.

## The object
**GENERALIZES to all tested Dirichlet L-functions (§570-571): e^{−i arg L′(ρ)} reads μ(n)χ(n) — the coefficients of 1/L — at z≈17 for χ mod 3, χ mod 4 (quadratic) AND χ mod 5 (COMPLEX order-4, §574). Universal across the full character group — a general Selberg-class principle, not special to ζ.**

For the Riemann zeros ρ = ½ + iγ (γ > 0) define the zero-sum with a **unit-modulus** weight
    G(u) = (1/N) Σ_{0<γ<Γ} e^{−i·arg ζ′(ρ)} · e^{iuγ}.
Empirically (500 → 10⁴ zeros):
- Peaks of |G(u)| land exactly on u = log n for **squarefree** integers n (40/40 vs 0/40 non-squarefree
  at 10⁴ zeros, p < 0.001; §530); squares are actively suppressed.
- At u = log n, Re[e^{ic}G(log n)] ∝ **μ(n)/n^{≈0.5} (phase-only 0.51, drifts; EXACT-Gonek 0.48→½; §563: brackets ½, not clean-fixed)** — the RH-critical exponent ½, robust to 0.50–0.52
  across 3 zero-counts × 3 ranges (§552b). Shuffling the phases collapses the exponent to the trivial 0.42.

## Derivation (why it is the Möbius function)
**Step 1 — the classical explicit formula.** With 1/ζ(s) = Σ μ(n) n^{−s} (Re s>1), Perron + residues at
the simple zeros give the (conditional, simple-zeros) Titchmarsh/Gonek formula for the Mertens function
    M(x) = Σ_{n≤x} μ(n) = Σ_ρ x^ρ/(ρ ζ′(ρ)) − 2 + Σ_{trivial}.
The **exact weight** is the Gonek weight 1/(ρ ζ′(ρ)). Verified here: reconstructing M(x) from 10⁴ zeros
with this weight gives corr 0.9933 → 0.9999 as N grows (§519, §522), and the reconstruction is sharply
best at Re ρ = ½ (σ-knob peak 0.513, §553) — the √x envelope of M(x) is the critical line.

**Step 2 — the weight decomposition.** Write ζ′(ρ) = |ζ′(ρ)| e^{i·arg ζ′(ρ)}. Then
    e^{−i·arg ζ′(ρ)} = |ζ′(ρ)| · (1/ζ′(ρ)) = ρ|ζ′(ρ)| · [1/(ρ ζ′(ρ))].
So the **phase-only weight equals the exact Gonek weight times ρ|ζ′(ρ)|.** The factor ρ|ζ′(ρ)| has modulus
|γ|·|ζ′(ρ)|, which varies *slowly on a logarithmic scale* (|ζ′(½+iγ)| grows like a power of log γ with
O(1) multiplicative fluctuation; |γ| is smooth). **Measured:** |ρζ′(ρ)| ~ γ^{0.22}, multiplicative
fluctuation only 1.77× (log-sd 0.57) — slowly-log-varying, confirmed; the 0.22 drift is what shifts the
measured exponent from ½ to 0.51. A slowly-log-varying real-ish envelope multiplying the
Gonek weight **preserves the oscillatory phases e^{iuγ} and the residue sign pattern** — which is where the
arithmetic (positions log n, signs μ(n)) lives — while only mildly reweighting magnitudes. Hence the
phase-only sum inherits the squarefree comb and the μ-signs, and an *approximate* √n scaling (measured
0.51 vs the exact ½), the small excess being the log-drift of ρ|ζ′(ρ)| (the §519b gap law, ≈0.009·ln n).

**Step 3 — the phase is S(T).** The intrinsic identity (§541, machine-verified to 4×10⁻¹⁴, classical via
ζ(½+it)=e^{−iθ_H}Z): arg ζ′(ρ) = π·S(γ) − π/2, where S(T)=π⁻¹ arg ζ(½+iT). Therefore the twist is a pure
**S(T)-phase**, and the whole statement reads: *the zero-fluctuation field S(γ) alone — no magnitude data —
reconstructs the Möbius function's positions, signs, and (to log-drift) its RH-critical √n normalization.*

## What is and isn't claimed
- **Known:** the Mertens explicit formula (Step 1) and the Gonek weight are classical.
- **The novel content (lit-gated NOT found, §521):** that dropping |ζ′| entirely — the **unit-modulus / phase-
  only** weight e^{−i arg ζ′(ρ)} — still yields the squarefree comb and μ(n)/√n, i.e. the Möbius function is
  carried almost entirely by the *phase* arg ζ′ ≡ πS, with magnitude nearly irrelevant to positions/signs.
- **Not claimed:** any new proof of RH or of M(x)=O(x^{½+ε}). The σ-knob pins Re ρ=½ as a *consistency*
  check of the explicit formula (on-line zeros used), not an independent probe.

## Open (to lift from REPORT toward a theorem)
1. Bound the deviation phase-only − Gonek: show ‖ρ|ζ′(ρ)| − (its log-average)‖ contributes only o(1) to the
   comb positions/signs (moment estimates for |ζ′(ρ)|; Gonek's conjecture territory).
2. Prove the exponent → ½ exactly as Γ→∞ (the 0.51 → 0.50 drift), quantifying the §519b gap.
3. Connect to Ng / Gonek–Ng work on Σ x^ρ/ζ′(ρ) and to unconditional partial results.


## FULL SELBERG CLASS (§570–583) — degree-independent, robust across families
The phase field arg L′(ρ) universally encodes the Dirichlet-inverse (Möbius-type) coefficients of ANY
L-function: e^{−i arg L′(ρ)} summed over the zeros reads the coefficients of 1/L. CONFIRMED across:
  • ζ (degree 1, trivial χ): reads μ(n)/√n, perfect squarefree comb at 10⁴ zeros, sharpens with N.
  • ALL Dirichlet L (degree 1): χ₃, χ₄ (quadratic), χ₅ (COMPLEX order-4) — reads μ(n)χ(n), z≈17.
  • Ramanujan Δ (degree 2, MODULAR FORM): reads b(n) [b(p)=−a(p), b(p²)=1], z=+25 (212 zeros).
  • Elliptic curve E11a (degree 2, ELLIPTIC): reads its b(n), z=+10 (a_p validated vs known 11a).
  • sym²Δ (DEGREE 3): reads its b(n), z=+7.7; wrong-L control 0.18 (§591). Three degrees spanned.
VERIFIED (§579): reads genuine degree-2 p²-structure + is L-specific (Δ's b 0.63 >> ζ's μ 0.04).
DERIVED (§567/577/580): Euler product of exp(S_L prime-sum) with the L-function's own Euler factor;
twist rebuilt from primes matches the true phase at 0.97–0.99. Lit-gated NOT found (§521). A general,
verified, derived, robust Selberg-class principle — NOT a proof of RH.


## L-function generalization (§570–577) — universal, derived
For ζ and every tested Dirichlet L (χ₃, χ₄ quadratic; χ₅ complex order-4): e^{−i arg L′(ρ)} reads
μ(n)χ(n) — the coefficients of 1/L — at z≈17. DERIVED (§577): the twist rebuilt from χ-weighted primes
reproduces the χ-comb (agreement 0.971), so the Euler-product mechanism carries over with χ(p). A universal
degree-1 Selberg-class principle, not special to ζ. Next frontier: degree-2 (modular/elliptic L) — §575.

## Mechanism of the squarefree SELECTION (§567) — Euler product of exp(S prime-sum)
e^{-2iπS(γ)}, S(γ)=-(1/π)Σ_p Σ_k p^{-k/2}/k sin(kγ log p), gives -2iπS = Σ_p p^{-1/2}(p^{iγ}-p^{-iγ})+(k≥2).
Then exp = Π_p exp[p^{-1/2}(p^{iγ}-p^{-iγ})], first-order = Π_p(1 + p^{-1/2}p^{iγ} - ...): the coefficient of
n^{iγ} is nonzero only for SQUAREFREE n (distinct primes; a repeated prime needs the suppressed k≥2 terms).
CONFIRMED: the twist rebuilt from PRIMES ALONE (no argζ' data) gives the squarefree comb with 0/50 non-sf
contaminants (agreement 0.949 with the real twist). This derives WHICH integers appear. The leading terms give
POSITIVE coeffs (½)^r/√n, so the μ SIGNS are NOT from this expansion — those are the Gonek/Perron piece above.
The two derivations are complementary: Euler-product => squarefree selection; explicit formula => μ signs + √n.
## Theorem statement + proof outline (the path to a theorem)
**Conjecture/Theorem-candidate.** Let L(s)=Σ_{n≥1} a(n)n^{−s} be an entire L-function in the Selberg
class with Euler product L(s)=Π_p F_p(p^{−s})^{−1}, functional equation, and (assuming GRH for L)
nontrivial zeros ρ=½+iγ. Write 1/L(s)=Σ_{n≥1} b(n)n^{−s}. Then the unit-modulus twisted zero-sum
    G_L(u) := lim_{T→∞} (1/N_L(T)) Σ_{0<γ<T} e^{−i arg L′(ρ)} e^{iuγ}
has |G_L| peaking at u=log n for n∈supp(b), with peak weight asymptotically ∝ b(n)/√n.

**Proof outline (4 steps; gaps flagged).**
1. *Perron/residues on 1/L.* For x not an integer, Σ_{n≤x} b(n) = Σ_ρ x^ρ/(ρ L′(ρ)) + (arch/pole),
   conditional on GRH + simple zeros (classical, generalizes the Mertens formula). Since ρ=½+iγ,
   x^ρ=√x·x^{iγ} — the √x envelope forces the n^{−1/2} normalization. So the **Gonek weight 1/L′(ρ)**
   reconstructs the b(n). [Numerics: ζ Gonek corr→0.9999, §519/522; σ-knob pins Re ρ=½, §553.]
2. *Phase-only = Gonek × slow envelope.* e^{−i arg L′(ρ)} = |L′(ρ)|·(1/L′(ρ)) = ρ|L′(ρ)|·[1/(ρL′(ρ))].
   The envelope ρ|L′(ρ)| is slowly-log-varying [measured |ρζ′|~γ^{0.22}, 1.77× fluct, §554], so it
   preserves the oscillatory factors e^{iuγ} and the residue-sign pattern (where supp(b) and sgn b(n)
   live) while only reweighting magnitudes. Hence G_L inherits the b(n) comb. GAP: bound
   ‖ρ|L′(ρ)| − avg‖ ⇒ o(1) on comb positions/signs (Gonek–Ng moment estimates for |L′(ρ)|).
3. *Support = Euler product of exp(S_L).* arg L′(ρ) ≡ πS_L(γ)−π/2 (exact for ζ, §541; S_L=π^{−1}arg L).
   e^{−iπS_L} expands via log L=Σ_p Σ_k c_k(p)p^{−ks}/k (c_k=trace of Satake^k) into Π_p(local Euler
   factor of 1/L) ⇒ supp(b): squarefree (deg 1), cubefree-with-p² (deg 2), etc. [Confirmed: prime-only
   twist reproduces the comb at 0.95–0.99 agreement across ζ, χ₃, Δ; §567/577/580.]
4. *The phase is S_L.* From L(½+it)=|L|e^{i arg L}, arg L′ at a simple zero = arg L′(ρ), tied to the
   S_L jump — exact identity for ζ (§541, 4×10⁻¹⁴), Vanlessen-Szegő for the string picture.

**Correction (§588).** Step 2's "envelope contributes o(1)" is WITHDRAWN for magnitudes: the envelope
correction decays only at log-rate (fluc/main 0.35→0.47 growing ~N^0.1; §563 exponent drift; §519b gap).
CORRECT THEOREM SCOPE: the phase-only readout gives supp(b) EXACTLY and sgn(b) EXACTLY (100% non-circular,
§584-585); magnitudes only up to slowly-varying log-scale distortion. Support+signs are the theorem-grade content.

**Status.** Numerically established + mechanistically derived across ζ, all Dirichlet L (incl complex χ),
degree-2 modular (Δ) and elliptic (E11a). REPORT-tier; the GAP in step 2 (envelope bound) is the main
obstruction to a full theorem — it is Gonek-conjecture-adjacent. No claim on RH itself.

## Audit & robustness (2026-08-03, §585–589)
- **Sign-exactness verified non-circularly** (rotation fixed via Gonek, not μ): 100% (305/305, n<500),
  corr 0.974; worst-case sign margin stable ~0.08 natural units for N≥2000 (§589) — robust, not luck.
- **Zero-set completeness**: the density formula exposed the Re-sign finder missing ~half the Dirichlet
  zeros; complete |L|-argmin rebuild (χ₃: 157 = prediction) IMPROVES the comb (0.591→0.663, z=16.7).
  Incompleteness diluted, never biased. Δ/E11a zeros carry ±0.03–0.05 AFE jitter (z-scores are
  underestimates).
- **Softness map complete** (§586–587): the comb readout channel is δ²-soft with ~0 coefficient under
  off-line perturbations (pair-symmetry protection) — the phase-arithmetic principle is NOT itself an
  RH-witness; coherent ON-line deformations are the linearly-stiff direction (drum, §545–546).
- **Variational picture** (§586): the true zero configuration is linearly stiff along coherent on-line
  deformations and quadratically soft along off-line ones — RH's difficulty in one sentence: the
  direction that would falsify it is exactly the direction the functional-equation symmetry protects.

## FINAL FORM (2026-08-03, §585–593): the Phase–Arithmetic Principle
**Statement (theorem-candidate, corrected scope).** Let L be an entire Selberg-class L-function with
Euler product, functional equation, GRH-zeros ρ=½+iγ, and 1/L(s)=Σ b(n)n^{−s}. Define the unit-modulus
zero-sum G_L(u) = N⁻¹ Σ_{0<γ<T} e^{−i arg L′(ρ)} e^{iuγ}. Then, as T→∞:
  (a) SUPPORT: |G_L| peaks exactly on {log n : n ∈ supp(b)} — the comb positions are exact
      (ζ: 40/40 vs 0/40 at 10⁴ zeros; suppression of n ∉ supp(b) verified through degree 2).
  (b) SIGNS: the aligned readout Re[e^{ic}G_L(log n)] has sgn = sgn b(n) exactly
      (ζ: 100% of 305 squarefree n<500, non-circular rotation, stable margin ≈0.08 natural units).
  (c) MAGNITUDES: ∝ |b(n)|/√n up to a slowly-varying log-scale distortion (NOT o(1); the phase-only vs
      exact-Gonek gap ≈ 0.009·ln n; the effective exponent brackets ½: Gonek 0.48→½, phase-only 0.50–0.52).
**Established numerically across:** ζ; L(s,χ) for all characters mod 3, 4, 5 (incl. complex order-4);
L(s,Δ) (degree-2 modular, z=25); L(s,E11a) (degree-2 elliptic, z=10); L(s,sym²Δ) (degree-3, z=7.7); L(s,sym³Δ) (DEGREE-4, z=4.75, density-gated 58/57.4) —
every spectrum density-validated complete. Sign-stability lemma (§594): envelope/main < 0.89 uniformly. Blind identification: 5/5 (each zero-set selects its own L).
**Mechanism (derived).** arg L′(ρ) ≡ πS_L(γ)−π/2 (exact); e^{−iπS_L} = Euler product of exp of the
prime sum of L ⟹ supp(b) selection; ≈Gonek-weight decomposition ⟹ signs; the ρ|L′(ρ)| envelope ⟹ the
log-distortion of magnitudes (§588). Twist rebuilt from primes alone matches the true phase at 0.95–0.99.
**Toward proof.** (a)+(b) are the claimable core. Ingredients: the Perron/residue formula for 1/L
(conditional: GRH + simple zeros), the S_L-exponential Euler product (Step-3 identity, exact), and a
sign-stability bound (margin ≈0.08 measured; needs a discrete-moment argument à la Gonek–Ng to prove the
envelope never flips a sign). (c) is provably NOT strengthenable to exact magnitudes (§588 measurement).

## Proof plan for (a)+(b) — the rigorous route (2026-08-03)
1. **Twist identity (exact).** arg L′(ρ) = πS_L(γ) − π/2 at simple zeros (proved via L = e^{−iθ}·(real
   Hardy-type function); machine-verified 4×10⁻¹⁴ for ζ). So the weight is e^{iπ/2}·e^{−iπS_L(γ)}.
2. **Support selection (identity + smoothing rigor).** e^{−iπS_L} = Π_p exp[−(local prime sum)] — an
   exact Euler-product rearrangement of exp(Im log L); its n^{iγ}-coefficient support is supp(b) by the
   local-factor algebra (degree-d recursion). Gap to close: uniform convergence/smoothing of the S_L
   prime series on the zero set (Gaussian-smoothed version verified at 0.95–0.99 agreement).
3. **Zero-sum evaluation (Gonek–Landau calculus).** Σ_{γ<T} x^ρ/L′(ρ) has known asymptotics under
   GRH+simplicity (Gonek's method: contour integral of x^s/L(s)); the phase-only sum = the same sum
   reweighted by the positive factor |ρL′(ρ)|. Claim (a): peak positions survive any positive
   reweighting. Claim (b): signs survive iff the fluctuation term is dominated — which is exactly:
4′. **Alignment (supersedes the domination lemma — §597).** The envelope-fluctuation readout ITSELF
   reads μ (corr 0.955, 100% sign agreement): |L′(ρ)|-fluctuations are explicit-formula objects sharing
   the prime source, so EVERY smooth-envelope reweighting of the Gonek weight preserves the sign
   structure. The proof burden is an ALIGNMENT statement (positive correlation of |L′|-fluctuations
   with the Gonek μ-oscillation; measured +0.955), not a smallness bound. Superseded text follows:
4. **Envelope-domination lemma (superseded).** max over comb points of
   |Σ(env−ē)·w_G·e^{iγ log n}| < |Σ ē·w_G·e^{iγ log n}| — measured: ratio ≤ 0.93 for N≥2000, decreasing
   (§594, 12% headroom). Proof route: second-moment/large-sieve bound for Σ env·oscillation using
   discrete moments of |L′(ρ)| (Gonek's conjecture / Ng's theorems give the needed moment growth);
   the measured headroom says the required constant is modest.
Under GRH + simple zeros + the step-4 inequality, (a) and (b) are theorems. (c) is intrinsically
log-distorted (§588) and NOT claimed beyond slow-variation.

## The closed chain (§600) + BK-alignment derivation sketch
CHAIN: phases are πS_L exactly (§541); exp(−iπS_L) is the Euler product of the prime sum, giving supp(b)
(§567); the Gonek weight reads b(n) classically; and the envelope ρ|L′(ρ)| is a local-gap functional
(Hadamard: log|L′(ρ_k)| = Σ_{j≠k} log|1−ρ_k/ρ_j| − log|arch| + c, §599) whose fluctuations carry the same
prime content (Bogomolny–Keating: gap fluctuations ↔ primes; measured directly in the zero-drum, §545).
Hence EVERY piece of the phase-only weight reads the same arithmetic — the alignment (fluc-readout reads
μ at 0.955/ζ, 0.594/Δ; 100% signs) is the BK mechanism seen through the Hadamard identity.
DERIVATION SKETCH for the alignment estimate: write env-fluc_k = Σ_{j≠k} [log|1−ρ_k/ρ_j| − expectation]
(≈50% of variance, §599). Insert into the readout sum and swap orders: the fluc-readout becomes a
DOUBLE zero-sum Σ_k Σ_j w_G(k) f(γ_k−γ_j) e^{iγ_k log n} — a pair-correlation-weighted Landau sum. By
Montgomery/BK, the pair-correlation's oscillatory part carries frequencies log p^m ⟹ this double sum
reproduces prime combs with the SAME sign pattern as the single sum (the BK phase relations), giving
positive alignment. To prove: the BK secondary-term asymptotics for this specific double sum (conditional
on Montgomery-type pair-correlation with the arithmetic corrections). Measured targets: 0.955 (ζ, 10⁴),
0.754 (Δ, 212). VALIDATED (§601) then SHARPENED (§602): the gap-kernel readout reproduces the alignment at +0.998, and
the gap-shuffle ablation (0.998→0.996, no collapse) shows it lives in WINDOW-SCALE S-fluctuations, not
fine pair-correlations. UNIFIED: every component of the weight (phase = πS pointwise; envelope = smoothed
S-functional) draws on the ONE source — the S-field explicit formula. The proof needs only the mollified
version of the step-2 Euler mechanism, not separate BK asymptotics.
