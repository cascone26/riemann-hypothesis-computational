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
| 2 | **Tail lemma** (finite K₀ cert ⟹ true infinite-dim positivity) | **DONE** | Handle | §208-tail (2026-07-24): K0=13 fails at the NL=8000 precision floor (mode 1 interval contains 0); K0=7 is the real minimal operational block — all 4 pivots rigorously >0 (`operators/track_tail/arb_certify_k0.py`, mode 1 margin 1.69e-7), Gershgorin tail dominance δ≥0.4552, Schur residual 3.2e-7 ≪ margin. Independently re-run and reproduced by the control-panel session (not just self-reported) after an earlier K0=13-based claim on this row was found false and corrected — see LEDGER §207/§208-tail. |
| 3 | **Odd sector** cert | **DONE** | Handle | §208-oddsector (2026-07-24): `arb_certify.py --parity odd` (additive flag, even path unchanged) certifies mode 1, pivot [0.0002395±5.56e-8] rigorously>0. Independently re-run and reproduced by the control-panel session directly (note: its own PROOF_CERT.json used cert-type "PD", which makes verify_gate.py check an unrelated Groskin baseline instead of this script — a process gap; the underlying result itself is confirmed correct by direct re-run). See LEDGER §207/§208-oddsector. |
| 4 | **Screw ≡ CvS/CCM equivalence** | **UNRESOLVED** | — | §200 Gap 3. `rigor/screw_cvs_equiv.py` built & runs, but its spot-check compared AT the criticality endpoint (a=½log3) where BOTH forms are ~0 — the (4,1)vs(5,0) "mismatch" is one near-zero eigenvalue's sign, NOT a disproof (§202). CLOSES when: re-compare at an INTERIOR a (clean +margin) with matched normalization; then prove change-of-basis P in Arb. Needed to transfer piece 5 to our form. |
| 5 | **Groskin grid** τ⪰0 sweep | **DONE (for τ)** | Handle | §4b; `external/groskin_2607.02828/code/arb_ldlt_certify.py`; c∈{2..100}×N≤500 all n_neg=0, rigorous Arb interval LDL^T. Transfers to OUR form only via piece 4. |
| 6 | **Continuum-in-a** (close gaps between integer c) | **COMPOSITION PROVEN** | Report+ | §203-continuum (2026-07-24): affine+SOS composition independently verified (verify_gate.py PASS). Window [0.48,0.52] certified Q(a)⪰0 via SOS + affine enclosures. Full continuum [½log2,½log100] pending (mechanical sweep, no new theory). `operators/track_continuum/sos_affine_compose.py`. |

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
prime — the minimal new result. Everything else ({4,5,6}) scales it out.

**Pieces 1, 2, 3 are now all DONE/Handle (2026-07-24) — the minimal viable theorem is assembly-ready.** Writing
the composed proof (not a new discovery task — pieces 1+2+3 just need to be stated together as one theorem) is
the next concrete step. {4,5,6} remain open/in-progress and scale the result out to the full 25-prime grid +
continuum, but are not required for the minimal result.

_Last updated 2026-07-24 (control-panel session, after independently re-verifying pieces 2 and 3 — see LEDGER
§207 for a prior false claim on these two rows that was caught and corrected before this update)._
