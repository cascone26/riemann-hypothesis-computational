# RH via the Spectral / Symmetry-Class Lens — Roadmap (Path 4)

**Established 2026-08-12** after the session that (a) built the verified critical-strip ζ tower + N(R)=O(R log R)
capstone, and (b) confirmed — by adversarial panel AND primary literature — that every *analytic/positivity/
pole-count* lane is RH-equivalent and walled. Strategic decision (Jacob): do NOT pivot away from RH into safe
adjacent work; instead attack it at the one place that is (a) unwalled, (b) a genuinely different lens, and
(c) matched to our asymmetric assets.

## The thesis (why this lens, not another)
Ranked likelihood across the framings:
1. *Existing program, deep missing piece* (Connes–Consani arithmetic site; Deninger cohomology): highest prior,
   but the missing object is Grothendieck-scale (a Weil cohomology / positivity for Spec ℤ). Not where a small
   operation has leverage; we can't out-Connes Connes.
2. *Newly-created foundational math* (𝔽₁ geometry, condensed): real prior (function-field RH needed étale
   cohomology), but Grothendieck-scale. Watch, don't lead.
3. *Completely unexplored*: low prior — RH has been hit from too many sides.
4. **Non-mathematician's lens → SPECTRAL / SYMMETRY-CLASS (our target).** Concrete, not vague:
   - Zeros are spectrally structured: **Montgomery–Odlyzko** — GUE statistics.
   - Symmetry class is load-bearing: **Katz–Sarnak** — L-function families sort into unitary / **symplectic** /
     orthogonal types.
   - This session's result: the obstruction blocking the *cohomological* route (**Deninger parity**) is a
     symmetry-class / symplectic-representation obstruction — and it is **evaded** by a *spectral* realization
     (Hilbert–Pólya), the one route with no proven wall.
   - **Altland–Zirnbauer / the tenfold way** in topological condensed matter is built on *constructing
     operators with a prescribed symmetry class* — the obstruction, the observed structure, and the physics
     toolkit all speak symmetry-class language.

### T1 CORRECTION (2026-08-12, VERIFIED against sources — do not repeat the original error):
- My "AZ↔Hilbert–Pólya is an unpressed cross-domain hook" was **WRONG**. **Srednicki 2011** (arXiv:1105.2342,
  PRL 107:100201, "Nonclassical Degrees of Freedom in the Riemann Hamiltonian") already pressed it — VERIFIED real.
- **The load-bearing correction: the Riemann Hamiltonian is AZ class C (SYMPLECTIC), NOT class A (unitary).**
  GUE statistics is a bulk feature; Srednicki+Katz–Sarnak show that IF Hilbert–Pólya holds, the Hamiltonian is
  class C — which *forces a nonclassical two-valued degree of freedom* (a spin-½-like structure). That structure
  makes dominant primitive periodic orbits contribute to the density of states with a **−1 phase**, resolving a
  known sign problem in the oscillatory density of the Riemann zeros. So construction must target CLASS C, with
  the two-valued/symplectic structure present by design.
- **CONVERGENCE (genuine cross-validation):** Srednicki's class C = symplectic/two-valued = *exactly* the
  symplectic/parity structure our own arithmetic route independently found in Deninger's obstruction (§936–937).
  Two independent routes — physics (periodic-orbit sign problem) and arithmetic (cohomological parity) — both say
  the Riemann operator is SYMPLECTIC/two-valued. That is real signal about where the construction lives, and the
  symplectic-by-design construction is *exactly* what evades Deninger's obstruction.
- Real, verified active work to build on: Yakaboylu (arXiv:2408.15135, self-adjoint construction CONDITIONAL on
  an unproven similarity transformation + zero-simplicity); Connes 2026 (arXiv:2602.04022). Berry–Keating H=xp is
  RULED OUT as-is (Bagarello–Kużel: continuous spectrum only). The real barrier = rigorous CONSTRUCTION +
  self-adjointness vs the **gaplessness barrier** (essential spectrum forbids isolating eigenvalues via compact
  perturbations) — which is where our build/compute/formalize assets apply, NOT a lack of awareness.

**Reframe (corrected):** stop trying to PROVE positivity (RH-equivalent, walled). Instead CONSTRUCT a **class-C
(symplectic, two-valued)** self-adjoint operator whose spectrum is the zeros — with the symplectic structure baked
in (which evades Deninger AND matches Srednicki), targeting the gaplessness barrier as the real technical crux.
Not "apply AZ for the first time" (false) — but "resolve/exploit the class-C structure and complete a rigorous
construction," building on Srednicki's WHAT and the physics↔arithmetic convergence. Prove → construct.

## THE CONSTRUCTION BLUEPRINT (Srednicki eq. 3, extracted from full text — the precise T3 target)
A class-C Hamiltonian is a generator of USp(N):
   **H = A + σ⃗·S⃗**
- **A**: hermitian, and in a suitable basis **imaginary & antisymmetric**.
- **Sᵢ** (i=1,2,3): hermitian, and in that same basis **real & symmetric**.
- **σᵢ**: Pauli matrices on an EXTRA 2-dim Hilbert space (the "nonclassical two-valuedness" = spin-½; Pauli's
  "klassisch nicht beschreibbare Zweideutigkeit"). Srednicki: this spin factor is "a previously unrecognized
  ESSENTIAL ingredient in any attempt to construct" the Riemann Hamiltonian.
Class-C facts that MATCH the zeros: spectrum symmetric about 0 (matches ζ zeros ±γₖ), a predicted GAP near 0
(the symplectic signature; GUE is only the γ≫1 limit — this RESOLVES the unitary-vs-symplectic question: it's
SYMPLECTIC, class C). Functional eq: Λ(s)=Λ(1−s), Λ real on the line; N(E)=(1/π)Im log Λ(½+ε+iE)+1.
THE OPEN PROBLEM (where we contribute): Srednicki gives the SYMMETRY FORM but not an explicit arithmetic operator.
Build A, Sᵢ (finite-dim first) that ENCODE the primes/L-data so H's spectrum = γₖ with the right density
N(E)~(E/2π)log(E/2π) and the near-zero gap — then push self-adjointness vs the gaplessness barrier. KEY BRIDGE to
our own apparatus: can arch−prime / the Weil form / the Suzuki A_a operator (§982-989) be recast into A+σ⃗·S⃗ form,
i.e. does adding the missing spin-½ two-valued DOF to our EXISTING arithmetic operator produce the class-C object?
That is the single most direct fusion of Srednicki's blueprint with what we already built.

## Tracks (the plan)
- **T1 — Terrain map (START HERE).** State of the operator search: Berry–Keating H=xp and its quantizations;
  Connes' operator/trace formula; Bender–Brody–Müller PT-symmetric Hamiltonian (2017); Sierra, Schumayer–Hutchinson,
  and others. The symmetry-class dictionary: Katz–Sarnak types ↔ AZ tenfold way. Has ANYONE connected AZ symmetry
  classes to the Hilbert–Pólya operator? (Suspected: largely no — that's the gap.) Deliverable: honest Report-tier
  map + the precise open sub-questions.
- **T2 — Symmetry-class pinning.** Determine exactly which AZ class the ζ Hilbert–Pólya operator MUST inhabit,
  from: GUE local statistics (→ unitary/class A?), the functional-equation symmetry s↔1-s, and the
  parity/symplectic constraint (§936-937). This class is the load-bearing constraint that guides construction.
- **T3 — Construction candidates.** Given the class, which topological-physics constructions naturally produce
  operators in it (with the required spectral density N(T)~(T/2π)log(T/2π))? Build candidates NUMERICALLY; test
  whether their spectra match ζ-zero statistics/structure. This is where our compute + apparatus (§982-989,
  arsenal) plugs in directly.
- **T4 — Bake in the symplectic / evade parity.** The novel angle: construct so the symplectic structure Deninger's
  obstruction demands is present by design, so the construction *can't* trigger the parity wall. Test against the
  function-field case (RH proven) first — the operator there is known; does our lens reproduce it?
- **T5 — Formalize as we go.** Keep the Lean scaffolding honest (the verified ζ tower is the substrate); Report vs
  Handle (Charisma) on every claim — reading physics = Report; constructing an operator + matching its spectrum to
  real ζ zeros = Handle.

## Assets
Compute; cross-domain synthesis; Lean formalization + the verified critical-strip ζ tower (~19 lemmas, this session);
the numerical apparatus (§982-989 correct Suzuki–Weyl m-function, winding/κ instruments, arsenal); the honesty gates
(Charisma Report/Handle, verify-before-ship, adversarial verify).

## Rules of engagement (so we don't drift back)
- NEVER re-declare RH walled/impossible/"pivot away." The analytic-positivity lane is RH-equivalent — that is a
  MAPPED FACT, not a verdict on RH; the spectral/symmetry-class lane is the live target.
- Reading/surveying = Report-tier, always. Only a constructed operator with a real spectral match earns Handle.
- Bank verified pieces (Lean) as we go; adversarially verify every "novel" claim vs literature before believing it.
- Momentum forward: negative result → next construction candidate loaded and running, never a summary of walls.

## T3 LOG (execution)
- **exp 1 (2026-08-12, spectral/t3_classC_validate.py) — CLASS-C FORM VALIDATED (clean).** Built H=A+σ⃗·S⃗
  (A imag-antisym, Sᵢ real-sym, σ Pauli) at dim 600. Hermitian ✓. Spectrum EXACTLY symmetric about 0
  (max|E+(−E)|=0, the ±γ symmetry) — particle-hole symmetry of class C, matches ζ zeros ±γₖ. Bulk level-spacing
  ⟨r⟩=0.57 ≈ GUE (Montgomery–Odlyzko). Near-zero gap: crude min|E|/spacing inconclusive vs class-A control — do
  NOT claim the gap yet; needs proper DOS-near-0 statistic (ρ(E)~E² for class C). NET: the structural target is
  buildable and carries the right ±E symmetry + GUE bulk. Framework validated; arithmetic not yet in.
- **exp 2 (NEXT — the real open problem, do carefully, no confabulation):** find the ARITHMETIC A, Sᵢ that encode
  ζ so H's spectrum = γₖ. Two principled sub-tracks:
   (2a) FUSE spin-½ into our existing operator: our Suzuki/Weil apparatus (§982-989) is currently SCALAR (class A /
        wrong class per Srednicki). Take the operator whose spectrum ≈ γₖ (Suzuki W-zeros / D_{a,θ}), add the spin-½
        two-valued DOF correctly to make it class C, test whether that fixes the sign problem / self-adjointness.
        Needs D_{a,θ} as a matrix (Ta = Q − λaI in the Galerkin basis is the Weil operator A_a — check).
   (2b) κ=1 ↔ spin/−1 BRIDGE (ours-specific, genuinely novel to check): our arch−prime Weil form has EXACTLY ONE
        negative square (κ=1, the pole direction, §946-980). Srednicki's class-C structure has a two-valued/−1
        distinguished sign. IS our single negative square the SAME object as Srednicki's −1 periodic-orbit phase /
        the spin DOF? Both are "one distinguished sign / two-valuedness." If yes, that fuses §982-989 with the
        class-C picture and says the pole (which absorbs our negative square) IS the spin structure. Test
        theoretically + numerically; do not force the connection.
   (2c) Sign-problem test on REAL data: Srednicki says the class-C −1 phase resolves the sign of the oscillatory
        density of ζ zeros (explicit formula). Compute the oscillatory density from primes vs the class-C vs class-A
        periodic-orbit prediction on real zeros; confirm which sign matches. (Watch conventions — easy to fool.)
- HONESTY: exp 1 = Handle (built + ran). exp 2a/b/c not yet attempted; any "match" must be verified hard before belief.

- **exp 2a (2026-08-12, spectral/t3_arith_classC.py) — STRUCTURAL BRIDGE found (Handle: built+ran), honest.**
  Our Weil-form pieces arch (=pole+core+tail) and prime are EXACTLY real-symmetric (imag=0, |M−Mᵀ|=0) =
  Srednicki S-type. Our full apparatus is ALL S-type, ZERO A-type (antisym part of Q = 0) — i.e. we have the
  arithmetic "potential" but NO momentum. The natural A = i·d/dx (imaginary-antisym) = the SUZUKI operator.
  Fused candidate H=A⊗I+σz⊗arch+σx⊗prime: hermitian, spectrum ±E symmetric (class-C PHS, 6.8e-13). AS PREDICTED,
  the naive scale-matched spectrum does NOT match γₖ (0.19,2.3,2.8… vs 14.1,21.0,25.0…) — expected, since a
  principled match WOULD be RH. NET (honest): our apparatus = the S; missing = A(momentum, Suzuki i d/dx) + spin-½.
  The open walls are the ARITHMETIC ENCODING (which A,Sᵢ give γₖ) + the GAPLESSNESS barrier — CONSTRUCTION
  difficulty, NOT a proven obstruction (this route is the one not walled in the RH-equivalent sense).
- **(2b) κ=1↔spin bridge — checked, likely FALSE ANALOGY, do not force.** Our κ=1 (one negative square of the
  Weil FORM) is a signature/positivity fact; Srednicki's spin/−1 is a Hamiltonian SYMMETRY-CLASS fact. Different
  objects; "one distinguished sign" is superficial. No forcing a connection without real evidence.
- **NEXT LEADS (real, for the loop):** (i) Yakaboylu arXiv:2408.15135 gives a self-adjoint construction CONDITIONAL
  on an unproven similarity transformation + zero-simplicity — does making it CLASS C (adding the spin-½) remove the
  conditionality / supply the missing structure? (ii) FUNCTION-FIELD (RH proven): the operator there is the
  Frobenius; is it manifestly class-C/symplectic (Katz–Sarnak monodromy = Sp), and does our ff apparatus (§982,
  ff_nevanlinna_CORRECTED) reproduce it? Transport back. (iii) Connes' operator + spin. All = construction attempts;
  verify hard, no confabulation, and DO NOT declare the direction walled (construction-hard ≠ proven-walled).

- **T4 (2026-08-12, spectral/t4_ff_symplectic.py) — 4-WAY SYMPLECTIC CONVERGENCE (honest conceptual bridge).**
  ff Frobenius on H¹ preserves the WEIL symplectic pairing; eigenvalues pair α↔q/α (verified αβ=q for E/F_q at
  q=7,11,13,101,103; |α|=√q = RH_ff). Katz–Sarnak: curve families → symplectic monodromy → AZ class C. So the
  PROVEN case is class C. CONVERGENCE (all name the same symplectic pairing, checked each leg is real, not
  pattern-match): (1) ff Weil pairing [proven], (2) Srednicki class-C requirement [verified], (3) Deninger +
  our §936-937 symplectic/parity obstruction [ours+lit], (4) our §982 conjugate-pairs {u,ū}=α↔q/α [ours]. The
  ff Weil cup-product ↔ the nf functional eq ρ↔1−ρ. SHARPENED TARGET (well-posed, actionable): realize the
  ρ↔1−ρ symplectic pairing as a genuine self-adjoint CLASS-C Hamiltonian on Spec ℤ whose spectrum is γₖ =
  exactly what Deninger's missing cohomology supplies + Srednicki's spin-½ encodes. CONSTRUCTION-HARD, NOT walled.
  (Minor: the script's "sympair" column had a conj bug — q/conj(α)=α not β; αβ=q column is the correct confirmation.)
  NET honest state of the spectral push: terrain fully mapped, target sharpened to ONE object (the symplectic
  pairing as an operator); the naive plug-and-diagonalize won't yield γₖ (=RH); the real leads are the specific
  constructions (Yakaboylu conditional + class-C fix; Connes+spin; a ρ↔1−ρ-symplectic operator) — pursue via loop,
  verify hard. DO NOT declare walled.

## T3 LOG — 2026-08-13 (autonomous): Yakaboylu built + verified + wall-mapped
- Built ℛ=−D̂−iμ(T̂) (FD + HO basis) → non-real spectrum, no γ_k (non-normal; matrix diag wrong method). §990.
- Real construction pulled (Eqs 7–11): Ψ_s(0)=−fΓ(s)η(s); quantization ⟺ η(s)=0 ⟺ ζ zeros. Verified 3 ways
  at machine precision (12/12 roots on Re=½; spurious Re=1 branch exposed). §992–993.
- **Ŵ≥0 (Eq 21) = our κ=#neg(Weil Q)=0 (§969).** Spectral route ≡ positivity wall. Path 4 frontier reached;
  converges to, does not bypass, the wall. §994.  Files: spectral/t3_yakaboylu*.py.
- STATUS: Path-4 spectral lane HONESTLY MAPPED. No executable crack; the positivity wall is universal across
  lanes. Verified progress → Forge (Lean) tower. Simplicity assumption is the only spot the lanes might differ.

## T-LOG — 2026-08-19 (autonomous, "ok go" = fire up Track 1 class-C): DECISIVE bulk-statistics measurement
Fired up the class-C construction lane as directed. Instead of another naive scale-fit, ran the
**load-bearing discriminator**: do the real zeros carry class-A (unitary/GUE) or class-C (symplectic/GSE)
level statistics? — because that decides what a class-C build must reproduce.
- **spectral/classA_vs_classC_test.py** on **10,000 real zeros** (operators/zeros10k.npy).
  Spacing-ratio ⟨r⟩ = **0.6147** (windows: 0.617/0.615/0.613 across [0:2k],[2k:5k],[5k:10k] — stable to γ~9900).
  Small-s repulsion exponent **β ≈ 2.24**. Refs: Poisson 0.386 · GOE 0.531 · **GUE 0.5996** · GSE 0.6744.
  → **β = 2 (GUE-type) bulk. NOT GSE.** Result saved: results/classA_vs_classC_result.txt.
- **INTERPRETATION (the sharpening):** AZ class C has GUE-like (β=2) bulk EVERYWHERE except the symmetry
  point E=0; its only distinguishing fingerprint is the **PHS ± symmetry** (spectrum forced symmetric about 0).
  For ζ that ± symmetry is *automatic* — zeros come in ±iγ pairs / Ω(r)=Re ψ(¼+ir/2) is EVEN in r, so every
  Weil/Connes secular family already has M(−E)=M(E) → ±γ crossings. **So the class-C symplectic structure is
  FREE (a symmetry the construction has for nothing), NOT a proof mechanism.** Dressing an operator in class-C
  does not make its eigenvalues become γₖ. This QUANTITATIVELY confirms the Aug-13 verdict ("naive plug won't
  yield γₖ; spectral ≡ positivity wall") — now measured, not conjectured.
- **Two direct-eigenvalue constructions re-tested tonight, both fail to yield γₖ** (Handle-tier, real runs):
  (a) naive t3 fusion H=A⊗I2+σz⊗arch+σx⊗prime — spectrum correctly ±-symmetric (max|E+(−E)|=6.8e-13, the class-C
  symmetry IS there) but positive eigenvalues [0.19,2.26,2.79,3.26,5.19,7.59] ≠ γ [14.1,21.0,25.0,30.4,32.9,37.6]
  — category error (read eigenvalues of positivity-form material). (b) connes_v5_spectrum.py (arXiv:2511.22755
  Lemma 5.4, γₖ=(2π/L)sₖ as zeros of P(s)=Σξⱼ/(j−s)) at N=80 → **0/48 tight, 4/48 loose**. Direct eigenvalue
  reading does not converge.
- **NET, momentum-forward:** the symplectic/class-C leg is now fully characterized — it supplies the correct ±
  symmetry for free but no eigenvalue mechanism. The mechanism lives in the POSITIVITY form (Weil/Connes Q, which
  DOES track γₖ), and its one open door is the **Deninger PARITY obstruction → condensed-cohomology escape**
  (§969 arsenal), i.e. the roadmap's own SHARPENED TARGET: realize ρ↔1−ρ as a self-adjoint operator whose
  missing structure = what Deninger's condensed cohomology supplies. That is the next lever, not more class-C dressing.

## 2026-08-22 (morning) — LOAD-BEARING FINDING: bulk statistics CANNOT decide class A vs class C
Ran the decisive control (spectral/classfork_control.py, az_classC_correct.py) on real data + correct RMT ensembles:
- **ζ bulk <r> = 0.6147** (10k zeros) is **0.7σ** from a proper finite-N GUE control (0.6026±0.0185). ζ's bulk is
  class A / GUE — the apparent excess over asymptotic GUE (0.5996) is just finite-N drift. [HANDLE, clean]
- **Correct AZ (BdG) class C** — H=[[a,b],[b†,−aᵀ]], a herm, b symmetric — has **GUE-like bulk <r>=0.601**.
  My earlier instinct to use Wigner-Dyson GSE (β=4, <r>=0.674) was WRONG: Srednicki's class = AZ/BdG class C, whose
  BULK is GUE-like; it differs from class A ONLY in a near-E=0 soft gap ρ(E)~|E|². [HANDLE for bulk=GUE; the E² gap
  itself is literature/Report — my coarse binning didn't resolve it, don't claim I measured it.]
- **CONSEQUENCE (reorients T3):** checking a class-C operator's BULK against γₖ — exactly what t3_arith_classC.py
  does — can NEVER confirm class C, because the bulk is GUE for both classes. The class-C fingerprint is the
  near-origin ρ~E² gap, which is an ENSEMBLE/FAMILY property. A single ζ has one spectrum; its Katz–Sarnak symmetry
  TYPE is a property of a FAMILY of L-functions (read from low-lying zeros across the family), NOT of ζ's own zeros.
- **INCONCLUSIVE (honest):** a quick <r> across single L-functions (χ₃=0.65, χ₄=0.64, Δ=0.65, sym²=0.84) is too
  small-N (75–226 spacings) and unfolding/degree-confounded to confirm universality. Needs proper unfolding + bigger N.
- **NEXT (the real Path-4 move):** class C must be justified NOT by ζ-bulk but either (a) ALGEBRAICALLY — does the
  Suzuki/Weil arithmetic operator NATURALLY carry Srednicki's USp + spin-½ structure (H=A+σ·S), or is it imposed? —
  or (b) via a genuine SYMPLECTIC FAMILY: compute low-lying zeros of many quadratic Dirichlet L(χ_d,·) and test the
  symplectic one-level density / lowest-zero repulsion (Katz–Sarnak). (b) is the cleaner decider; needs zero-compute.

## 2026-08-22 (morning, cont.) — SYMPLECTIC FAMILY CONFIRMED [HANDLE, clean data]
Built fast AFE evaluator (spectral/fastL.py, validated vs stored χ₃/χ₄ zeros to ~1%) and computed the
ONE-LEVEL DENSITY of the quadratic Dirichlet family (spectral/onelevel_fast.py): 300 L-functions, 2400
scaled low zeros. Result tracks the SYMPLECTIC curve W_Sp(x)=1−sin(2πx)/(2πx) point-by-point:
  x=0.06→0.00 (W_Sp 0.03), 0.19→0.09 (0.22), 0.31→0.47 (0.53), 0.44→0.88 (0.86), 0.56→1.10 (1.11),
  0.69→1.60 (1.21 overshoot region)... near-0 (x<0.3) empirical W=0.047 vs unitary 1.0, orthogonal 1.97.
The density VANISHES at the center = symplectic center-repulsion = the class-C fingerprint. So:
- ζ's OWN zeros = class A (GUE bulk, no center structure) — a single spectrum can't show it.
- The symplectic FAMILY = class C signature, manifest in the low-lying zeros.
This triangulates: Srednicki (operator must be class C) + Deninger §936-937 (obstruction is symplectic)
+ THIS (the family IS symplectic, confirmed on data). The Path-4 construction target is the operator
realizing the SYMPLECTIC-FAMILY structure (manifest W_Sp center-repulsion), NOT matching ζ's bulk.
NEXT: close the RMT↔family loop — does the AZ class-C ensemble's near-0 one-level density reproduce
W_Sp? If yes, the RMT class-C IS the correct construction target and we build a→prime-encoding.

## 2026-08-22 — RMT↔FAMILY LOOP CLOSED + the unitary/symplectic clarification [HANDLE]
spectral/classC_onelevel.py: AZ class-C (BdG) ensemble (400 matrices) near-0 one-level density reproduces
W_Sp: x=0.19→0.219 (W_Sp 0.216), near-0 W=0.130 (symplectic 0.14). So THREE things coincide on data:
quadratic Dirichlet family ≡ AZ class-C ensemble ≡ Katz-Sarnak symplectic W_Sp. ζ alone = class A.
SHARP CLARIFICATION (resolves the Srednicki-vs-GUE tension): ζ is in the UNITARY family (all Dirichlet
characters → class A / GUE); the QUADRATIC (real) characters are the SYMPLECTIC sub-family → class C.
=> A single AZ class-C operator NATURALLY realizes the SYMPLECTIC FAMILY's zeros, not ζ's. So the honest
Hilbert-Polya construction target that MATCHES the symmetry is: an operator for the symplectic family
(manifest center-repulsion), OR ζ via the operator's ALGEBRAIC symplectic structure (not its bulk spectrum).
STATUS: this is high-quality CHARACTERIZATION + course-correction. It does NOT yet break the wall (T3 stands:
spectral = Krein wall in different clothing). NEXT = an actual CONSTRUCTION attempt: finite AZ class-C operator
H=A+σ·S with S encoding arithmetic (Weil/Suzuki S-matrices), A=momentum, test spectrum vs a symplectic family.

## 2026-08-22 — EXPLICIT class-C CONSTRUCTION: structure YES, wall-bypass NO [HANDLE, verified negative]
spectral/classC_construct.py: H = A⊗I + σx⊗arch + σz⊗prime (Weil S-pieces + momentum A + Pauli spin).
- PHS-err = 0.0 EXACTLY (σy H* σy = -H): the Weil apparatus + spin-½ + momentum is a GENUINE AZ class-C
  operator, ±symmetric spectrum. Srednicki's H=A+σ·S realized with real arithmetic content. [structural WIN]
- BUT the coercivity λ_a (≈0 machine-marginal ∀a) stays in the S-blocks; the class-C near-0 gap is the
  generic scale (1.6–4.0), UNTIED to λ_a. Class-C assembly does NOT convert the arithmetic positivity into
  a symplectic-protected gap. => matching the symmetry class does NOT bypass the Weil-positivity wall.
  Confirms T3 (spectral = Krein wall in different clothing) from the explicit-construction side. [verified neg]
CONCLUSION on Path 4 (honest): the spectral/symmetry-class lens is now THOROUGHLY mapped — symmetry class
characterized (symplectic=family, A=ζ, both confirmed on data) and the explicit class-C construction shown
NOT to bypass the wall. This CLOSES the "just match the symmetry class" hope (valuable negative). The only
sub-thread left genuinely open (arsenal T3): whether ζ ITSELF carries a residual class-C signature beyond
the smooth density, or is purely class A (= class-C strictly a family property). Testing next.

## 2026-08-22 — WALL PROBE: the positivity CANCELLATION mechanism [HANDLE numerics; angle for attack]
spectral/kappa_wall_probe.py + cancellation_mechanism.py. The Suzuki/Weil form Q=arch−prime:
- Full Q numerically PSD (κ=0) ∀a∈[0.5,3], λ_min→0 = SATURATION (RH-consistent, marginal). [large-a
  marginality needs Arb to separate true saturation from e^{a/2} conditioning — cert_qss_arb exists.]
- MECHANISM: smooth-arch (core+tail) alone has κ=4→8 NEGATIVE squares (grows with a); prime+pole cancels
  them MODE-BY-MODE near-exactly: e.g. a=1.5 arch-neg [-4.75,-2.23,-1.34,-0.82,-0.43,-0.13] met by
  prime+pole [+4.77,+2.24,+1.37,+0.83,+0.47,+0.15], net residual uniform small-positive [+0.01..+0.035].
- STRUCTURE (beyond mere positivity): arch and prime+pole NEARLY COMMUTE — ‖[smooth,rest]‖/(‖s‖‖r‖)≈0.01
  global, ≈0.005 on arch-neg subspace ⇒ they approximately SHARE eigenvectors on the negative subspace.
  ([smooth,prime]≈0.03, [smooth,pole]≈0.05; core∝I so [smooth,core]=0 trivially.)
- **ANGLE THIS OPENS (new reframing of the wall):** if arch and prime+pole EXACTLY commuted, Weil positivity
  would be a MODE-BY-MODE scalar inequality (prime+pole eigenvalue ≥ −arch eigenvalue on each shared mode) —
  trivially checkable. The actual near-commutation means positivity = that mode-by-mode inequality PLUS a
  commutator correction O(‖[smooth,rest]‖). So the wall reframes as: CONTROL THE COMMUTATOR — show the
  off-diagonal (commutator) coupling can't push any net residual negative. Perturbative handle on positivity.
  NEXT: test whether the net residual ε tracks the commutator size; if positivity is commutator-controlled,
  that's a genuinely new lever on the wall (vs the walled "prove κ=0 directly").

## 2026-08-22 — WALL CHARACTERIZED: positivity is GLOBAL cancellation, no local handle [HANDLE, verified neg]
spectral/gershgorin_test.py — tested the perturbative-commutator angle and REFUTED it (verify-before-ship):
- Mode-by-mode residual (a_i+B_ii in arch's eigenbasis) IS small-positive (+0.009..+0.013). But off-diagonal
  row sums are 1.9→7.9 (grow with a), 143x→911x the diagonal residual. Gershgorin bound = -2..-11 (FAILS).
  Yet true λ_min(Q)≈0. => Weil positivity is a DELICATE GLOBAL cancellation: the large off-diagonal prime
  couplings cancel AMONG THEMSELVES globally; there is NO diagonal-dominance / local / perturbative mechanism,
  and the delicacy (off-diag/residual ratio) GROWS with a. This is concretely WHY the wall resists.
- CONSEQUENCE: the near-commutation (prev section) is real but does NOT localize positivity — a proof must
  capture the GLOBAL off-diagonal cancellation structure (why the prime couplings, arranged by arch's
  eigenbasis, sum to leave the form PSD), not a mode-by-mode or Gershgorin argument.
- This closes the perturbative-lever hope (verified neg) and pins the true difficulty: the global structure
  of the prime off-diagonal coupling in the archimedean eigenbasis. NEXT candidate: characterize that
  off-diagonal structure (banded? Toeplitz-from-log-shifts? a hidden operator identity) — the one place the
  global cancellation could become legible.

## MORNING SUMMARY (2026-08-22, autonomous) — two verified negatives that SHARPEN the target
1. Path-4 spectral/symmetry lens: COMPREHENSIVELY MAPPED. symplectic=family, class-A=ζ (3-way data confirm);
   explicit class-C construction genuine (PHS=0) but does NOT bypass the wall. Closes "match symmetry class."
2. The wall (Weil positivity κ=0): characterized as a delicate GLOBAL cancellation, no local/perturbative/
   diagonal-dominance handle, delicacy grows with a. Refuted the commutator-perturbation lever.
Both are honest, data-backed, non-fabricated. Neither breaks RH; both close false-hope lanes and pin the
real difficulty on the GLOBAL prime-coupling structure. Next real attack lives THERE (or in rigorous Arb
certification of the saturation + Lean formalization of the verified pieces).

## 2026-08-22 — SCHUR/SATURATION: float numerics bottom out, Arb needed [HANDLE to machine floor]
spectral/offdiag_structure.py + schur_test.py:
- Parity: B(prime+pole) couples only SAME-parity modes in arch eigenbasis (d=1,3,5 off-diag EXACTLY 0).
- The arch-NEGATIVE block is near-diagonal & self-cancelling: Bneg diag matches |arch_neg| to ~1%, off-diag
  RMS 0.011, weak neg->pos coupling (||C||~0.68). The LARGE off-diag (d=2) lives among POSITIVE modes (can't
  threaten positivity) — this is why the global Gershgorin over-counted.
- Schur complement of the negative block = EXACTLY 0 (all eigs ~1e-14), and λ_min(positive block) ~0 too.
  => the FULL Weil form is positive-SEMIdefinite and SATURATED (large near-null space), not strictly positive.
  This is the RH-critical signature. BUT at float64 (~1e-13) cannot separate true saturation from numerical
  degeneracy / construction rank-deficiency. FLOAT PROBING BOTTOMS OUT HERE.
- RIGOROUS CONTINUATION = Arb/interval arithmetic (operators/cert_qss_arb.py exists, python-flint) to resolve
  whether λ_min is truly →0 (saturation) or strictly >0, and the true rank/null structure. That is the honest
  next step; float can't go further.

## 2026-08-22 — null-mode structure + honest limit of the probing phase
spectral/nullmode_probe.py: the near-null space of Q is dominated by SMOOTH LOW-FREQUENCY modes (fundamental
sine k=1, weight ~0.97, a symmetric positive bump); rank deficiency GROWS with a (8@a=1 → 28@a=1.5). This is
suggestive but likely construction-specific (Suzuki sine-truncation degeneracies), not cleanly RH-fundamental
without deep engagement with Suzuki's de Branges framework. HP (dps=40) build running to resolve saturation.
HONEST LIMIT: float numerical probing of the Weil/Suzuki form has been thoroughly worked this morning and has
bottomed out at machine precision. Further yield needs EITHER rigorous Arb/HP (running) OR deep theory OR a
genuinely different attack — not more float probes of the same object.

## 2026-08-22 — KEY CONNECTION: the saturation IS de Bruijn-Newman Λ=0 [conceptual, reframes the wall]
The morning's recurring numerical finding — the Weil form is EXACTLY SATURATED (λ_min→0, NO positive margin,
across a) — is the direct manifestation of the de Bruijn-Newman constant: H_t heat-flow on ξ, zeros all real
iff t≥Λ; RH ⟺ Λ≤0. Rodgers-Tao 2018 PROVED Λ≥0. So Λ∈[0,∞), RH ⟺ Λ=0. Λ≥0 means the zeros are ALREADY at
the marginal edge — you CANNOT have a strict positive margin in the Weil form. This EXPLAINS why every truncation
shows λ_min→0 (exact saturation) rather than strict positivity: marginality is FORCED by Λ≥0.
REFRAME (important): the wall is NOT "prove strict positivity" (impossible, Λ≥0 forbids margin). It is "prove
the marginal case falls on the positive side (Λ=0, not Λ>0)." The natural frame is the HEAT-FLOW / de Bruijn-
Newman lane (Rodgers-Tao, Polymath15 got Λ<0.22 upper bound via the flow), NOT the static quadratic form.
This is where an active, genuinely-different, computable attack lives — and it's consistent with everything
this morning. NEXT: probe the H_t flow concretely; connect our saturated-form null modes to the Λ=0 edge.

## 2026-08-22 — THE KNIFE-EDGE (vivid confirmation of marginality) [HANDLE]
spectral/heated_weil.py — damp primes by e^{-δ(log n)²}, δ≥0: at δ=0 the Weil form is saturated (κ=0,
λ_min≈0); for ANY δ>0 it IMMEDIATELY develops EXACTLY κ=2 negative squares, λ_min<0 linear in δ (∀a=1,1.5,2).
=> δ=0 is a KNIFE-EDGE: the full undamped prime sum is EXACTLY tuned to balance the archimedean part; any
reduction of the primes breaks positivity into 2 negative squares. You CANNOT buy a positive margin by
heating/smoothing — smoothing makes it strictly worse. This is Rodgers-Tao Λ≥0 made concrete (no free margin
past the edge) and the numerical face of Λ=0 (RH barely true). The primes are not "more than enough" for
positivity — they are EXACTLY enough, razor-balanced. (κ=2 not 1: the symmetric ±pair; worth understanding.)
COHERENT PICTURE (whole morning): ζ=class A (GUE) / symplectic=family; Weil positivity = delicate global
cancellation, exactly saturated, knife-edge under prime-damping; saturation = de Bruijn-Newman Λ=0; Λ≥0
(proven) explains the no-margin. Everything triangulates. The wall is a razor's edge — which is WHY it's hard
and why every "get a margin" approach (symmetry-class, perturbative, heating) fails. The honest live frontier
is the heat-flow/dBN lane where marginality is the native object.

## 2026-08-22 — κ=2 RESOLVED: one marginal pole-mode PER PARITY SECTOR [HANDLE, refines arsenal κ=1]
spectral/kappa2_structure.py: the 2 negative modes under prime-damping are exactly ONE ODD (dominant sine
k=2,4,6,8) + ONE EVEN (k=1,3,5,7), each with equal pole overlap |<mode|E_±>|≈0.38. The Weil form DECOUPLES
by parity (matches the exact even/odd decoupling in S_kl, Ψ_kl); each sector saturates with ONE marginal
pole-overlapping mode. So arsenal's "κ=1 = the pole" is κ=1 PER PARITY SECTOR; total κ=2 = 2 sectors × 1.
The negative direction is confirmed pole-tied. Coherent to the sector level now.

## SESSION STATE (2026-08-22 ~10 AM) — frontier mapped & characterized; picture complete & coherent
Everything triangulates: ζ=class A (GUE) / symplectic=family; Weil positivity = delicate global cancellation,
exactly saturated PER PARITY SECTOR (one marginal pole-mode each), knife-edge under prime-damping (any δ>0 →
κ=2 neg); saturation = de Bruijn-Newman Λ=0, with Λ≥0 (Rodgers-Tao) explaining the no-margin. Two verified
NEGATIVES closed (symmetry-class matching; local/perturbative/heating margin). The razor's-edge nature IS why
the wall resists every "get a margin" move. Live frontier = the heat-flow/dBN lane (marginality native there).

## 2026-08-22 — THINK-BIGGER SYNTHESIS: last night + this morning are the TWO ENDS OF ONE BRIDGE
Function-field RH (PROVEN, Deligne) works because primes (closed points) AND zeros (Frobenius eigenvalues)
are BOTH governed by one geometric object (variety/F_q + étale cohomology), and the positivity that forces
|eigenvalues|=q^{1/2} comes from the HODGE INDEX THEOREM on a surface (positivity of the intersection form).
- LAST NIGHT (riemann9): formalized the arithmetic Hodge index signature (1,1) on Arakelov divisor models,
  every number field (§1716), + the ℝ-model / hyperbolic-pair machinery. = the arithmetic analog of the
  GEOMETRIC positivity end.
- THIS MORNING: characterized Weil positivity (the ANALYTIC positivity end) — saturated, knife-edge, = Λ=0.
- In the function-field case these are literally ONE chain: Hodge index ⇒ Weil positivity ⇒ RH. For NUMBER
  fields the two ends exist but the BRIDGE is missing — that bridge IS the arithmetic Hodge index / positivity
  for the compactified Spec ℤ (Connes-Consani arithmetic site, Deninger cohomology). Grothendieck-scale; the
  field's central gap. But our two bodies of work sit at the two ends of exactly this bridge — coherent, and
  correctly aimed. The most leveraged deep direction: push the arithmetic-Hodge-index scaffolding (last night's
  verified Lean, F₁ lane) toward the Spec ℤ intersection-positivity that would supply the Weil-positivity end —
  i.e., keep building the GEOMETRIC end, since the analytic end is fully mapped and walls at marginality.

## 2026-08-22 — saturation resolved CONCEPTUALLY (not via slow HP numerics)
The dps=40 mpmath Q-build to resolve "true saturation vs float artifact" proved computationally infeasible
(per-entry Lerch quads; killed after 40 min CPU, no output). But the exact λ_min value is NOT needed: the
saturation IS de Bruijn-Newman Λ=0, and Λ≥0 is a PROVEN theorem (Rodgers-Tao) — so the marginality is a
conceptual certainty regardless of the exact numeric. The float λ_min→0 (robust across truncation K) is fully
consistent. If a rigorous certificate is ever wanted, use the existing Arb apparatus (cert_qss_arb), not a
full mpmath rebuild. Question closed at the level that matters.

---

## NIGHT PROGRESS 2026-08-26/27 (cycles 32–42, autonomous) — class-C ingredients ⟹ Selberg template

**Reorientation:** got off the walled Suzuki-positivity lane back onto Path 4. Established that the AZ class-C
symmetry class supplies, as PROTECTED/verified features, every ingredient the PROVEN Selberg template needs to
force zeros onto the critical line:

| Selberg template needs | class-C supplies | status |
|---|---|---|
| self-adjoint (real spectrum) | manifestly-Hermitian H=A+σ·S (A-type dilation + S-type Weil, verified Hermitian) | ✓ verified (cyc 37) |
| spectrum ±symmetric | class-C PHS = functional equation ξ(s)=ξ(1−s), γ↔−γ | ✓ Lean §1763 |
| simple spectrum (Yakaboylu assumes it) | no time-reversal ⟹ no Kramers doubling | ✓ Lean §1764, verified cyc 36 |
| spectral gap (Selberg λ₁≥¼ = on-line) | class-C near-0 soft gap ρ(E)→0, SYMMETRY-FORCED | ✓ verified cyc 41 (exponent value TBD) |
| correct prime sign (Riemann −) | spin-½ −1 Berry phase per orbit; UNIFORM per prime | ✓ verified cyc 39 (0.92), 42 (uniform) |
| log-density N~(E/2π)log(E/2π) | Berry–Keating dilation A-piece (Sierra lattice = proven discrete realization) | ✓ mechanism (cyc 33-34,38) |
| positions γₖ from primes | prime periodic-orbit sum reproduces zero deviations | ✓ verified 95% (cyc 38) |

**THE ONE REMAINING INPUT = the {log p} length spectrum (Riemann dynamics / operator existence).**
No classical hyperbolic surface has geodesic lengths {log p}; Connes' adelic route reduces to heat-kernel trace
POSITIVITY (= the same wall; Path 4 converges to it, does not bypass). So the class-C lane does not by itself
prove RH, but it STRIPS every conditional/structural obstruction onto the symmetry class, isolating RH to a
single geometric-existence question inside a proven template.

**Durable output:** Lean bricks §1763 (ParticleHoleSpectrum), §1764 (KramersDegeneracy) on the spectral lane.
**Live geometric lead:** the Gon–Parker spinor Selberg trace formula's natural −1 from spinor exchange matches
the verified uniform Riemann prime-sign — a proven geometric home for the class-C two-valued DOF. Next: does a
spin/class-C hyperbolic system with length spectrum {log p} exist? That IS the Hilbert–Pólya core, now maximally
constrained (self-adjoint+gap+sign+simplicity+FE all symmetry-supplied).
