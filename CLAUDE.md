# Riemann Project — Agent Context Rules

## Subagent briefing — lean by default (toggle: say "full context" to override)

PROCESS.md (196KB) and STATUS.md (129KB) are huge and MUST NEVER be trimmed or
overwritten — `feedback_process_provenance.md` requires the full cumulative
history, append-only. But every subagent (Task/Agent tool call) is its own
isolated context that can't share the main session's cache — if it reads the
full PROCESS.md+STATUS.md to orient itself, that ~150-200K-token read is paid
fresh, from scratch, every single time. Launching 3 agents in parallel to
fetch 3 papers, found 2026-07-15, cost ~3x that fixed tax on top of the
actual (much smaller) paper-fetching work.

**Default rule:** when launching a subagent for a narrow, independent task
(fetch/summarize a specific paper, check a specific claim, run a specific
calculation) — do NOT tell it to read PROCESS.md or STATUS.md. Give it only
the specific inputs the task needs (the paper URL/topic, the specific
question, the specific claim to check). The main session already has the
full history cached — it does the synthesis after the subagent reports back.

**Toggle — when a subagent genuinely needs full history:** if the task
requires deep prior context (e.g. "review our whole approach and find a
contradiction," "continue the proof from where we left off"), explicitly say
so and it gets the full PROCESS.md/STATUS.md briefing. That's the exception,
not the default — ask for it by name rather than getting it automatically.

**Why this matters:** same research pace, same depth, same number of
parallel agents when genuinely needed — this only removes the *duplicate*
cost of each one independently re-reading the entire process history from
scratch when its actual task didn't need it.

## Multi-session protocol (2026-07-23) — applies whenever more than one
session/agent is working this repo concurrently

Read `SESSIONS.md` and `ATTACK_QUEUE.md` at the start of every session — not
just `LEDGER.md`. If you are one of several concurrent tracks (you'll know
because you were briefed with a track name and a worktree path), the
following is mandatory, not optional:

1. **Claim before you start.** Your item should already be marked
   `<-- CLAIMED by:<track> since:<ts>` in `ATTACK_QUEUE.md`. If you're
   popping a new item (yours is exhausted), claim it there first, and note
   it in `SESSIONS.md`'s roster.
2. **Stay in your lane.** Only create/edit files inside your own
   `operators/track_<name>/` and `results/track_<name>/`. The shared
   append-only docs (`LEDGER.md`, `PROCESS.md`, `STATUS.md`,
   `ATTACK_QUEUE.md`, `SESSIONS.md`) — always append at the very end, never
   edit or reformat an existing entry, even a stale one. Editing an existing
   *shared code file* (e.g. `arb_certify.py`, `rigor/rigor_lint.py`,
   `engine/verify_gate.py`) is allowed but must be logged in `SESSIONS.md`'s
   shared-file table BEFORE you edit, and must be strictly additive (new
   function/flag/mode — never rewrite existing behavior). If your edit would
   change existing lines another track might also need, stop and treat it as
   a conflict (see 5), don't just proceed.
3. **Sync often, not at the end.** Run `./sync-track.sh <name> "<msg>"` from
   your worktree after every real finding (every LEDGER append) — not batched
   at session end. This is the actual mechanism by which you see what the
   OTHER tracks have found (the sync script rebases `main` back into your
   branch), so syncing often is how you avoid duplicating work, not just
   hygiene.
4. **No unverified breakthrough claims — the Charisma Report/Handle gate is
   mandatory here.** A numeric result is REPORT. It only becomes HANDLE
   (CERTIFIED / PROVEN language, an ntfy breakthrough ping, marking an
   `ATTACK_QUEUE.md` item done-by-proof) after `engine/verify_gate.py`
   independently re-verifies it — run it and cite its PASS/FAIL output in the
   LEDGER entry. This already exists in this repo; the only new rule is that
   it's non-optional in the multi-session setting, where an unverified claim
   from one track can misdirect the other two.
5. **Real conflicts go to Jacob, not to a guess.** If `sync-track.sh` reports
   a conflict, it already aborted the merge/rebase and logged it to
   `SESSIONS.md`'s NEEDS JACOB section — don't try to resolve it yourself;
   stop that line of work and flag it plainly when you next report.
6. **Self-drive across your lane.** When your claimed item is exhausted,
   check `SESSIONS.md` + `ATTACK_QUEUE.md` for the next unclaimed relevant
   item before stopping — same anti-stall discipline as always
   (`rh-working-rules`: never end on "what next?", keep the next method
   already loaded, no stopping-off-ramp language, ETAs when asked).
7. **Flag stalling, don't grind silently.** If you've tried the same
   sub-approach ~5 times without a genuinely new result (not "5 runs
   confirming the same thing" — 5 real *attempts* at cracking it), update
   your `SESSIONS.md` roster status to `stalled: <N> attempts on <what>, no
   new result` instead of quietly continuing or switching approaches
   unannounced. This is not a request to stop — keep working — it's a
   signal for a Fable path-check (Jacob's PRISM brainstorm review), which is
   triggered by this kind of friction, not run continuously.
