# SESSIONS — multi-track coordination hub

Live roster for concurrent RH sessions. Read this (and `ATTACK_QUEUE.md`) before claiming work. Full
protocol lives in `CLAUDE.md` ("Multi-session protocol"). One canonical checkout stays on `main` at
`~/projects/riemann` — never edit files there directly while tracks are running; it's the integration
point only.

## Roster

| track | claimed item | worktree | branch | status | last sync |
|---|---|---|---|---|---|
| tail | Tail lemma finalize (§201 -> K0 cert) | ~/projects/riemann-tracks/tail | track/tail | starting | - |
| continuum | a-interval feasibility (§199.6/199.7) | ~/projects/riemann-tracks/continuum | track/continuum | starting | - |
| oddsector | Odd sector cert + screw/CvS equivalence check (§200 gap2/gap3) | ~/projects/riemann-tracks/oddsector | track/oddsector | starting | - |

Update the `status`/`last sync` cells (or append a new roster row block below if the table gets unwieldy)
each time `sync-track.sh` runs. Status values: `starting`, `running`, `blocked`, `done-item` (popped next
queue item), `stalled` (see NEEDS JACOB).

## Shared-file touch log
Append here BEFORE editing any file outside your own `operators/track_<name>/` or `results/track_<name>/`
namespace (e.g. `arb_certify.py`, `rigor/rigor_lint.py`, `engine/verify_gate.py`). Edits to shared files
must be strictly additive (new function/flag/mode — never rewrite existing behavior).

| when | track | file | why |
|---|---|---|---|

## NEEDS JACOB
Empty unless a real merge conflict or ambiguous call comes up. A track that hits a real conflict aborts
the merge (does not guess a resolution) and appends here with: track, files, both commit SHAs, one-line
description of the actual disagreement.

(none yet)

## Sync log
Append-only. `sync-track.sh` writes one line per sync attempt here — this is how a track sees what the
others have done without polling. Don't edit past lines.

(none yet)

## Protocol summary (full version in CLAUDE.md)
1. Claim your item in `ATTACK_QUEUE.md` (`<-- CLAIMED by:<track> since:<ts>`) before starting.
2. Work only in `operators/track_<name>/`, `results/track_<name>/`, and append-only additions to
   `LEDGER.md`/`PROCESS.md`/`STATUS.md`/`ATTACK_QUEUE.md`/`SESSIONS.md`.
3. Run `./sync-track.sh <name>` after every real finding (every LEDGER append), not just at session end.
4. No CERTIFIED/PROVEN/HANDLE language or ntfy breakthrough ping without running `engine/verify_gate.py`
   first and citing its PASS/FAIL output in the LEDGER entry.
5. When your claimed item is exhausted, check this file + `ATTACK_QUEUE.md` for the next unclaimed
   relevant item before stopping. Self-drive — never end on "what next?".
