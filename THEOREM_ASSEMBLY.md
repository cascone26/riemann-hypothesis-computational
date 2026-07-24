# THEOREM ASSEMBLY — live status of the reframed first theorem (§200)

**Read this + LEDGER.md at session start.** This is the narrow composition tracker: which pieces of the target
theorem are PROVEN, which are NUMERIC-ONLY, which are OPEN. Each row carries a proof pointer (Charisma
Report/Handle + the hardcoded proof-gate: "done" needs a proof pointer or it's written as claimed). Update the
status column when a piece moves; never mark a row Handle/PROVEN without a pointer to a re-runnable artifact.

## Target theorem (reframed)
> Weil positivity in **true (infinite-dimensional)** form, **unconditional**, for support a ≤ ½ log c for every
> integer c ≤ 100 (a ≈ 2.30, first 25 primes) — assembled from the tail lemma + the already-certified Groskin grid;
> then the continuum-in-a closes the gaps between integer c. Honest scope: finite-a. RH = the a→∞ limit, untouched.

## Pieces & status

| # | Piece | Status | Tier | Proof pointer / what closes it |
|---|-------|--------|------|-------------------------------|
| 1 | Finite even-block Arb cert at a point | **DONE** | Handle | §199.5; `operators/dbn_dynamics/arb_certify.py`; `results/weil_p2_cert_a0549.prov.json` (lint PASS, flint 0.9.0). All 4 LDL pivots rigorously >0. |
| 2 | **Tail lemma** (finite K₀ cert ⟹ true infinite-dim positivity) | **IN PROGRESS (de-risked)** | Report+ | §201/§202: exponents locked (diag ≈0.31·log kβ; coupling C_k=2k[kβI_k−g(2a)], g(0)=0). Stress-test found NO crack (tail λmin≥1.15, `rigor/tail_stress.py`); K₀≈13 sized, near-null Schur residual ~1e-7 (`rigor/k0_solver.py`). CLOSES when: Gershgorin δ(K₀) + Schur residual written as a proof with explicit constants, then arb_certify the K₀=13 block (modes [1,3,5,7,9,11,13]) at NL≥8000. |
| 3 | **Odd sector** cert | **OPEN** | — | §200 Gap 2. Easy (odd sector margin 4–40× wider than even) but arb_certify does even modes [5,7,3,1] only. CLOSES when: add even-k (odd-sector) modes to the certifier and re-run. |
| 4 | **Screw ≡ CvS/CCM equivalence** | **UNRESOLVED** | — | §200 Gap 3. `rigor/screw_cvs_equiv.py` built & runs, but its spot-check compared AT the criticality endpoint (a=½log3) where BOTH forms are ~0 — the (4,1)vs(5,0) "mismatch" is one near-zero eigenvalue's sign, NOT a disproof (§202). CLOSES when: re-compare at an INTERIOR a (clean +margin) with matched normalization; then prove change-of-basis P in Arb. Needed to transfer piece 5 to our form. |
| 5 | **Groskin grid** τ⪰0 sweep | **DONE (for τ)** | Handle | §4b; `external/groskin_2607.02828/code/arb_ldlt_certify.py`; c∈{2..100}×N≤500 all n_neg=0, rigorous Arb interval LDL^T. Transfers to OUR form only via piece 4. |
| 6 | **Continuum-in-a** (close gaps between integer c) | **IN PROGRESS** | Report | §199.6/§199.7/§202. `rigor/affine_arith.py` — affine primitive PROVEN (rigorously encloses a real entry, tighter than naive, advantage∝1/h). `rigor/sos_continuum.py` — SOS interval-cert MACHINERY PROVEN on a conditioned family. Composition identified: SOS must consume affine's sub-margin enclosure (a polyfit is ~1000× too coarse for the ~1e-6 margin). CLOSES when: affine applied to the full entry matrix → high-precision SDP certifies each window. Retires the panel-farming scheduler. |

## Composition logic (how the pieces make the theorem)
- **Point positivity, true form, at integer c:** piece 1 (finite cert, both parities = 1+3) **+** piece 2 (tail
  lemma) ⟹ true infinite-dim positivity at that a. Piece 4 lets piece 5's 99 certified c-points substitute for
  re-running piece 1 at each c.
- **Whole windows (continuum):** piece 6 (affine+SOS) upgrades each integer-c point to the open interval around it
  (between integer c the prime set is constant — smooth territory).
- **Assembled:** {1,2,3} give the machine's own proof-of-life at points; {4,5} extend it across 99 conductors for
  free; {6} fills the continuum. All finite-a.

## What "done" means for the FIRST theorem (minimal viable)
Pieces **2 (tail lemma) + 3 (odd sector)** applied to piece **1** at the window a∈(½log2, ½log3) = the p=2-only
positivity theorem, true form, both parities, unconditional. That alone extends Yoshida/Bombieri past the first
prime — the minimal new result. Everything else ({4,5,6}) scales it out. Current blocker for the minimal version:
finish piece 2 (tail lemma constants + K₀) and do piece 3 (odd sector — mechanical).

_Last updated 2026-07-23. Pieces 2,4,6 have tooling building under `rigor/`._
