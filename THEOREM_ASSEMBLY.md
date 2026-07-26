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
| 2 | **Tail lemma** (finite K₀ cert ⟹ true infinite-dim positivity) | **RECOVERED at K0=9 (§227); Report+/Handle-pending-tail-rigor** | Report+ | §208's K0=7 was WRONG (at K0=7 the near-null Schur residual 1.78e-7 > low-block margin 1.65e-7 — composition FAILS; §224/§226 caught this). But at NL=24000 (NL=1200 was a precision artifact) the tail COMPOSES at **K0=9**: low λmin 8.6e-8 > near-null Schur residual 4.8e-8, tail Gershgorin δ=0.86>0 (k0_solver, float). LOW-BLOCK now Arb-rigorous PD (9-mode superset cert, all pivots >0, NL=50000, §227). REMAINING for full Handle: Arb-interval version of the tail Gershgorin margin + Schur residual (k0_solver is float). See LEDGER §224/§226/§227. |
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

> **⚠ CORRECTION (2026-07-24 roundtable, Hawk) — piece 3 is VACUOUS and the line above is retracted.**
> `arb_certify.py --parity odd` certifies `modes=[1]`, but per §192 the ODD parity sector = EVEN-k sines
> (EVEN-v = odd-k). k=1 is an odd index — the cert certifies a 1×1 minor of the EVEN block. Receipt: its
> pivot 2.395e-4 = A·Q(1,1) = 0.549·4.3633e-4 exactly (the tool omits /A). The odd sector had NO certified
> entries. See `operators/track_tail/ODD_SECTOR_GAP.md`: first even-k measurements (Q(2,2)=0.0216, corner
> NOT diagonal-dominant), full transfer of the Loewner/u_c=−1/diagonal-law machinery verified at 1e-13, and
> a corrected cert (`arb_certify_oddparity_k0.py`, modes {2..24}, NL=50000) fired — result in
> `results/arb_certify_oddparity_K24_nl50000.log`. Piece 2's tail architecture (4′ cut) must ALSO be
> mirrored for the odd sector; same constants u_c=−1, β=4w_p/π, parity-flipped 1/m terms.

> **⚠ CORRECTION (§226, 2026-07-24 pm) — piece 2 (tail lemma) status DOWNGRADED to UNVERIFIED.** Re-running the
> verifying tool rigor/k0_solver.py shows it does NOT support "K0=7 ⟹ true-infinite-dim": it runs at NL=1200 (far
> too low — reports negative low-block λmin where arb_certify_k0 at NL=8000 shows PD), its own verdict is "No K0≤41
> satisfied both," and at K0=7 the near-null Schur residual (3.24e-7) EXCEEDS the low-block margin (1.69e-7) — the
> §208 "≪" was backwards. So the FINITE-BLOCK certs (pieces 1, 3 — K0=7 even block + odd, pivots >0 at NL=8000) are
> solid, but the INFINITE-DIM UPGRADE (piece 2) did not reliably compose at K0=7. **UPDATE (§227) — RESOLVED, the
> more-positive answer:** NL=1200 WAS the artifact; at NL=24000 the tail COMPOSES at **K0=9** (low λmin 8.6e-8 >
> near-null Schur residual 4.8e-8, tail Gershgorin δ=0.86>0). §208 simply used the wrong K0 — K0=7 genuinely fails
> (residual > margin), K0=9 works. The K0=9 low-block is now Arb-rigorous PD (§227). **Minimal theorem status:
> true-infinite-dim RECOVERED at K0=9 (Report+; full Handle pending an Arb-INTERVAL tail-bound — k0_solver is float).**
> The §208 "K0=7 DONE" record is corrected to K0=9; the theorem itself stands.

_Last updated 2026-07-24 (control-panel session, after independently re-verifying pieces 2 and 3 — see LEDGER
§207 for a prior false claim on these two rows that was caught and corrected before this update)._

## p={2,3} EXTENSION status (2026-07-24, LEDGER §217-218) — first certified crossing past the first prime
- **Multiprime cert tool built + gate-verified** (`operators/track_crossblock/arb_certify_multiprime.py`):
  generalizes the p=2-only screw cert to all prime powers q with log q<2a; reproduces the canonical p=2 pivots
  EXACTLY at a=0.549. The reusable industrialization tool.
- **CERTIFIED (Handle, Arb, VERIFIED §223/§225):** even-4×4 + odd sector both PD across a∈[0.549, 0.64]. Even sector
  certified at 12 points a∈[0.550..0.64], margin ~2-4e-8 (NLERCH=50000); odd sector certified across the full window
  [0.55..0.80] (§222). Above a=0.64 the even margin dips below ~1e-8 approaching the p=5 criticality (a=½log5=0.805),
  undetermined at NL50000 (§225). So both-parity reach = [0.549, 0.64] — first certified crossing past the first
  prime-power boundary, gated by the even margin's approach to the next criticality.
- **CORRECTION (§223):** an earlier reading (§217b/§218) that the even margin "erodes to sub-floor / fails past 0.57"
  was a NLERCH-precision artifact — the pivot RADIUS (~1/NLERCH² tail bound) was swamping the true margin. At
  NLERCH=50000 the margin resolves as a stable positive; a=0.62 (which "failed" at NL16000) certifies. Retracted.
- **HONEST SCOPE:** still FINITE-GALERKIN (even-4×4 + odd blocks), NOT yet true-infinite-dim — the p={2,3} mode-TAIL
  lemma remains open/delicate (margin ~3-4e-8 is small vs the §208-style tail Schur residual). So this is weaker than
  the minimal theorem (p=2, which HAS the tail §208). Symbolic closed forms for the entries: NOT found (§223) — the
  extension is by rigorous Arb certification, not symbolic proof.
- **RELAY route CLOSED (§218):** the norm-budget relay (the last operator-side/PLAN_PRIME_INDUCTION hope) is
  quantified — arch=−prime to 6 digits on the near-null (the prime budget tracks arch deepening near-exactly) —
  but it delivers only the LEADING cancellation; the actual margin is the RESIDUAL (~1e-8), = the prolate
  uncertainty defect (Fuchs e^{−4πλ²}), whose SIGN is the RH content. So the relay reduces to the residual/prolate
  wall like every operator-side route (§182). The operator-side is now multiply-confirmed exhausted; the forward
  path is the genuinely-new sign-carrying mechanism (§182 spec vi), not more operator-side certification.
