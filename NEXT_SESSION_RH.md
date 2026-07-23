# NEXT SESSION (updated 2026-07-23, LEDGER through §181)

## CURRENT STATE (§181) — every basin is mapped + walled; STOP basin-hopping (verified, not assumed)
Read the LEDGER, not this doc's older sections, for basin status — they are STALE. Verified closures:
operator-positivity (§122, six faces) · zero-density estimates (§180, route-3, not-an-RH-path) · mollifier proportion
(route-7, θ→1 = RH-hard) · **DBN Λ≤0 (§35-47e — a Λ-DETECTOR, no prime-handle; NOT "dormant/unmapped")** · **LeClair
spectral-flow (§125 — = Lagarias' known RH-equivalent)** · Jensen/moment/de Branges/Beurling-Nyman (OBSTRUCTION_MAP 1-5).
Do NOT re-enter these as if fresh (my §180 nearly re-treaded DBN; search-before-build/read §35-47e,§125 caught it).

**The honest forward set (§181), goal-relevant only:**
1. **§126 Lens B (the one concrete live sub-problem):** can Guth-Maynard's 2024 large-value/energy method (singular
   values of the prime matrix M_W) transfer to push the Levinson-Conrey mollifier past θ=4/7 (proportion past ~42%)?
   §180 verified GM's RESULT is pure density (doesn't touch the mollifier); the METHOD-transfer is open/possibly-
   unclaimed. Hard analysis, not a probe. Real progress if it works (proportion, not full RH unless θ→1).
2. **Non-classical cross-basin ideation** (§125 rule-1; superPRISM on Fable) — the only route to the genuinely-new
   idea the verified wall demands.
3. **Consolidate** the §176 cascade (novel analytical near-null anatomy, §175-flagged unpublished) + the maps.
Per riemann-goal (prove RH, not partials): 1 & 2 are the goal-relevant moves; 3 is banking.

---
_(older handoff below — LEDGER-through-§123, RETAINED for provenance but basin-status is SUPERSEDED by §181 above)_

# NEXT SESSION (updated 2026-07-21, LEDGER through §123)

## STRATEGIC PIVOT (§123): operator-positivity basin EXHAUSTED → enter the ZERO-DENSITY basin
Forks A+B (parallel) settled the strategy. FORK A (§122): attacking the limit-positivity directly is the RH-hard core,
no non-circular handle (arch O(−e^a) ↔ prime +O(e^a) cancellation is vivid but positive⟺reality⟺RH). Our operator-
positivity basin (§1-122, six equivalent faces §121) is exhausted at ONE wall. FORK B (§123, survey): the genuinely-
DIFFERENT basins (not our wall) are, ranked: ① ZERO-DENSITY/Levinson (obstruction=mollifier Cauchy-Schwarz losses,
GENUINELY different; ≥40.88% of zeros proven on-line — real unconditional progress; ACTIVE: Guth-Maynard 2405.20552
2024); ② Spectral-flow/LeClair 2406.01828 (S-matrix unitarity→real zeros, novel, unproven); ③ Hilbert-Pólya (density
mismatch); ④ dBN Λ≤0 (dormant). SAME-WALL (avoid): function-field/Connes arithmetic-site, moments/RMT, Beurling-Nyman.

## Immediate next block — enter the zero-density basin
1. **READ Guth-Maynard 2405.20552** (2024, new Dirichlet-polynomial large-value estimates → improved zero-density
   N(σ,T)) + the Levinson-Conrey mollifier method (how 41% on-line is proven; where the Cauchy-Schwarz loss is).
2. Find the concrete tractable sub-problem: improved Dirichlet-poly bounds, or a mollifier-free gap argument toward a
   higher critical-line proportion. This is a genuinely DIFFERENT obstruction (analytic, not positivity-at-limit).
3. High-risk alternative to keep in view: LeClair spectral-flow (rigorously justify S-matrix unitarity as Re(s)→½).

## What the operator-positivity campaign leaves (durable)
The comprehensive six-face map (LEDGER §121), tools (compact G_a §110, mpmath hi-prec §109, prime-bend §112, semilocal
Jacobi §119), and the honest closure: no operator-positivity face is a shortcut past the a→∞/ω→0/all-primes limit.
Don't re-enter this basin. Self-corrections banked (§92, §117): check definitional properties (|m|≡1) before overlays.

## Orientation
LEDGER §92-123, PROCESS/STATUS 07-21, memory rh-suzuki-screw-function-frontier. Sources in membrane_pulls/. Tools:
operators/dbn_dynamics/ (many). Model: Opus for the zero-density reads + work; Fable = brainstorm only.
Parked: REBOUND Λ-hunt; Rick offline (SpaceEngine specimen).

---
_(superseded detail below from the §106 redirect — kept for provenance)_

# NEXT SESSION (updated 2026-07-21 step-back, LEDGER through §117)

## State after the §117 retraction (read this first)
The §106 redirect (off the k_λ wall onto Suzuki's W/Cor-1.6 door) STANDS, but §113-116's "reality transfers for free
via Herglotz" was RETRACTED in §117: m_a = v̂_+/v̂_- is UNIMODULAR/INNER on R (|v̂_+|=|v̂_-| exactly, forced by Thm 1.5),
not Herglotz; inner limits do NOT forbid C+ zeros. The correct target in this framework is the classical one:
**E_ξ = Ξ + iΞ′ is Hermite–Biehler ⟺ RH** (equivalently Θ_ξ=E_ξ^♯/E_ξ inner in C+).

## The honest map: ONE wall, four faces (all equivalent, all = RH)
1. k_λ ≈ θ_x rigorously (CCM 2511.22755 §7 — the worldwide stuck point)
2. θ(a) selects the right self-adjoint extension in the a→∞ limit (Suzuki Cor 1.6; deficiency (1,1) = limit-circle)
3. E_ξ Hermite–Biehler / Θ_ξ inner (classical)
4. momentum→zeta spectral transition (the prime-bend §112 made qualitative progress on this face)
Everything this session touched lands on one of these faces. Do NOT re-derive the equivalences — they're mapped.

## What genuinely stands from §104-117
- Suzuki Thm 1.5: finite-a reality + simplicity UNCONDITIONAL (his theorem, source-read). The door is real.
- Compact G_a reformulation (§110, g_a_compact.py): breaks the large-a conditioning wall, clean λ_a>0 to a=3 float64.
- High-precision run (§109): §42's negativity = float64 artifact; even-simple robust; two floors separated
  (precision fixed; basis-resolution structural).
- Prime-bend (§112): primes inject zeta-like gap-irregularity (CV 0.393 vs zeta 0.389 at a=2), density-independent.
- S_t(z)/P_t(z) explicit (§110, from [14] eq 1.5-1.6) — the Parseval identity ∫S_x S̄_y dz = g̃(x,y) (RH-equivalent,
  §107b) is fully computable and NEVER attacked analytically. Genuinely untouched face.

## Next block (recommended order)
1. **READ CM 2310.18423** ("Zeta zeros and prolate wave operators" — semilocal prolate = arch+primes; queued twice,
   never read; the paper working the SAME arch↔prime structure as our §104).
2. **READ Suzuki 1204.1827** ("A canonical system of differential equations arising from ζ") — the explicit canonical-
   system Hamiltonian H(t) for Ξ; plus the de Branges subspace CHAIN (total ordering, 2606.09096 §5.7) — the chain IS
   the canonical system in inverse spectral theory, and we have NEVER explored it. HB(E_ξ) ⟺ H(t)≥0 is the classical
   reformulation — find where the primes sit in H(t).
3. Then: attack the **Parseval/completeness identity** (face nobody's worked): ∫S_x(z)S̄_y(z)dz = g̃(x,y) with S_t
   explicit — completeness/frame-theory tools, not approximation.
Method rule sharpened (§117 lesson): before claiming object X has property P, run P's DEFINITIONAL one-line check
first (|m|≡1 would have killed the Herglotz claim instantly). Magnitude checks before interpretive overlays.

## Orientation
LEDGER §92-117 (the full arc incl. retraction §117), PROCESS 07-20/21, memory rh-suzuki-screw-function-frontier.
Sources in repo: suzuki_2606.09096.txt (Thm 1.5 l.385, §7.2 l.2206, §8 l.2420), suzuki_2301.00421 (S_t eq 1.5-1.6),
rt_1801/, connes_2602/. Tools: operators/dbn_dynamics/ (validated: g_a_compact, g_a_W_prime_bend, mpmath_weil_hiprec;
confounded-documented: lp_lc_test, weyl_mfunction). Unpulled: 2310.18423, 1204.1827, 2112.05500.
Parked: REBOUND Λ-hunt (Olympus); Rick offline (SpaceEngine specimen unchecked).
Model: Opus for the reads + analytic work; Fable = brainstorm/step-back only.
