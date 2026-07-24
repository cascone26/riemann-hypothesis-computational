# SESSIONS — multi-track coordination hub

Live roster for concurrent RH sessions. Read this (and `ATTACK_QUEUE.md`) before claiming work. Full
protocol lives in `CLAUDE.md` ("Multi-session protocol"). One canonical checkout stays on `main` at
`~/projects/riemann` — never edit files there directly while tracks are running; it's the integration
point only.

## Roster

| track | claimed item | worktree | branch | status | last sync |
|---|---|---|---|---|---|
| tail | Tail lemma finalize (§201 -> K0 cert) | ~/projects/riemann-tracks/tail | track/tail | done-item | 2026-07-24T01:57:21Z |
| continuum | a-interval feasibility (§199.6/199.7) | ~/projects/riemann-tracks/continuum | track/continuum | done-item | 2026-07-24T01:58:17Z |
| oddsector | Odd sector cert + screw/CvS equivalence check (§200 gap2/gap3) | ~/projects/riemann-tracks/oddsector | track/oddsector | done-item | 2026-07-24T01:58:15Z |

Update the `status`/`last sync` cells (or append a new roster row block below if the table gets unwieldy)
each time `sync-track.sh` runs. Status values: `starting`, `running`, `blocked`, `done-item` (popped next
queue item), `stalled` (see NEEDS JACOB).

## Shared-file touch log
Append here BEFORE editing any file outside your own `operators/track_<name>/` or `results/track_<name>/`
namespace (e.g. `arb_certify.py`, `rigor/rigor_lint.py`, `engine/verify_gate.py`). Edits to shared files
must be strictly additive (new function/flag/mode — never rewrite existing behavior).

| when | track | file | why |
|---|---|---|---|
| 2026-07-24T02:35:00Z | oddsector | arb_certify.py | additive: add --parity arg (even/odd), select modes by parity, keep even path intact. Per protocol, logging before edit. |

## NEEDS JACOB
Empty unless a real merge conflict or ambiguous call comes up. A track that hits a real conflict aborts
the merge (does not guess a resolution) and appends here with: track, files, both commit SHAs, one-line
description of the actual disagreement.

(none yet)

## Sync log
Append-only. `sync-track.sh` writes one line per sync attempt here — this is how a track sees what the
others have done without polling. Don't edit past lines.

- 2026-07-24T02:40:00Z main (canonical) ok commit:9938b3f oddsector-piece3-verified (LEDGER §208 verify_gate PASS)

## Protocol summary (full version in CLAUDE.md)
1. Claim your item in `ATTACK_QUEUE.md` (`<-- CLAIMED by:<track> since:<ts>`) before starting.
2. Work only in `operators/track_<name>/`, `results/track_<name>/`, and append-only additions to
   `LEDGER.md`/`PROCESS.md`/`STATUS.md`/`ATTACK_QUEUE.md`/`SESSIONS.md`.
3. Run `./sync-track.sh <name>` after every real finding (every LEDGER append), not just at session end.
4. No CERTIFIED/PROVEN/HANDLE language or ntfy breakthrough ping without running `engine/verify_gate.py`
   first and citing its PASS/FAIL output in the LEDGER entry.
5. When your claimed item is exhausted, check this file + `ATTACK_QUEUE.md` for the next unclaimed
   relevant item before stopping. Self-drive — never end on "what next?".

### CONFLICT 2026-07-24T01:27:06Z — track:tail
merge track/tail -> main failed. Log: /tmp/sync-track-merge-tail.log

- 2026-07-24T01:27:06Z track:tail FAILED commit:cf3193b merge-conflict-see-tmp-log
- 2026-07-24T01:28:12Z track:tail ok commit:60a247a push-failed(non-fatal)
- 2026-07-24T01:29:52Z track:tail ok commit:6f741b6 push-failed(non-fatal)
- 2026-07-24T01:57:21Z track:tail ok commit:a25c05a push-failed(non-fatal)
- 2026-07-24T01:58:15Z track:oddsector ok commit:a6dbec4 push-failed(non-fatal)
- 2026-07-24T01:58:17Z track:continuum ok commit:aa25f36 push-failed(non-fatal)
- 2026-07-24T02:00:08Z track:tail ok commit:b6b30ab push-failed(non-fatal)

### CONFLICT 2026-07-24T02:10:02Z — track:tail
merge track/tail -> main failed. Your work IS still committed locally at 66bad18. Log: /tmp/sync-track-merge-tail.log

- 2026-07-24T02:10:02Z track:tail FAILED commit:66bad18 merge-conflict-see-tmp-log
- RESOLVED 2026-07-24T02:41:00Z by control-panel session: real conflict (both tail and oddsector landed
  concurrently, overlapping edits to LEDGER/ATTACK_QUEUE/SESSIONS/results/PROOF_CERT.json + an add/add on
  operators/track_tail/arb_certify_k0.py against an already-merged earlier version of the same path).
  Resolved by keeping both sides' real content (K0=7 fix + odd-sector cert), not discarding either.
  tail correctly stopped and reported per protocol rather than resolving it itself — this is what rule 5
  is for.
