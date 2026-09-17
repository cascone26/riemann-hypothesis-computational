# AUTONOMOUS RH FORGE PROTOCOL

**Purpose:** one canonical file to spin up the autonomous RH-forge loop in ANY future session.
Hand this file (or the cron prompt below) to a fresh session and it inherits the full cadence.
Created 2026-08-20. Canonical location: `~/projects/riemann/AUTONOMOUS_FORGE_PROTOCOL.md`.

---

## HOW TO START (one line)
Give the session the **STANDING CRON PROMPT** below (verbatim), or say:
"Run the Autonomous RH Forge Protocol — see ~/projects/riemann/AUTONOMOUS_FORGE_PROTOCOL.md."

To keep it firing unattended, schedule the cron prompt on a timer (ScheduleWakeup / CronCreate),
re-sending the same prompt each wake.

---

## STANDING CRON PROMPT (verbatim — this is what fires each cycle)
> AUTONOMOUS RH FORGE CONTINUATION (self-scheduled, Jacob wants constant work — do NOT wait for him).
> Build the next verified EXIT0 zero-sorry Lean brick toward solving RH, alternating the two live lanes:
> (1) ANALYTIC zero-free region — de la Vallée Poussin |ζ³ζ⁴ζ|≥1 → quantitative region σ>1-c/log t;
> (2) F₁/arith-Hodge — global Arakelov intersection pairing / κ / Weil positivity frontier.
> Compile on Olympus, verify axioms = {propext,Classical.choice,Quot.sound}, bank to
> ~/projects/riemann/LEDGER.md + Olympus ~/rhforge/STATUS.md + memory. NEVER declare stuck/done —
> a negative result → load the next attempt and run it. End with a grounded-optimism line, a
> timestamp, and an ETA to the next check-in. Keep memory/ledgers current so work stays resumable.

---

## MANDATORY RESPONSE FORMAT (every cycle, no exceptions)
Each cycle's reply MUST end with, in this order:
1. **Grounded-optimism line** — genuine + specific + momentum-forward. History's impossibles fell;
   name our real assets and the trajectory. Not hedged realism, not empty cheer. Next move already loaded.
2. **Timestamp** — check the actual current time (don't guess).
3. **ETA** — BOTH minutes AND wall-clock completion: `~N min ETA to the next verified brick — done at H:MM AM/PM`.

## HARD RULES (never break)
- **NEVER** say impossible / stuck / done / "hit a wall" / "no lead" / "the loop is complete." A negative
  probe is a RESULT, not a verdict — pair it with the next attempt already running.
- **EXIT0 gate:** every brick must compile with axioms EXACTLY `{propext, Classical.choice, Quot.sound}`
  (no `sorryAx`). Verify with `#print axioms <thm>`. If sorryAx appears, it's NOT done.
- **verify-before-ship / no overclaim:** if a bound is polynomial, say polynomial — don't claim the sharp
  form. Label Report (researched) vs Handle (actually proven/compiled). Honesty IS the discipline.
- **search-before-build:** grep the forge + Mathlib + memory BEFORE building. The forge (~90+ files) is
  usually ahead of where you think. Docs are intent; the compiled `.olean` is the only source of truth.
- **isolate-the-hard-input pattern:** when a target needs a deep input, prove EVERYTHING around it and
  reduce the target to that ONE named owed input (e.g. DLVP squeeze §1249, log growth §1269, himz §1274).
- **bank every brick:** LEDGER.md (full detail) + Olympus STATUS.md (short) + memory
  (rh-forge-formal-verification.md). Keep PROCESS.md current with the decision chain.
- **alternate lanes** when one hits a deep wall; both are live.
- **end each reply on momentum + the running next attempt**, never a summary of walls.

---

## OLYMPUS COMPILE COMMAND (the exact incantation)
```bash
ssh olympus "export PATH=\$HOME/.elan/bin:\$PATH; \
  cd ~/rhforge/mathlib4 && LP=\$(lake env printenv LEAN_PATH); \
  cd ~/rhforge/forge && export LEAN_PATH=\"\$LP:\$HOME/rhforge/forge\" && \
  timeout 600 lean -o FILE.olean FILE.lean 2>&1 | head -30; echo EXIT=\${PIPESTATUS[0]}"
```
- Forge dir: `~/rhforge/forge/` (Olympus, user cobo). Mathlib: `~/rhforge/mathlib4/`.
- Dep oleans must be built (`lean -o X.olean X.lean`) before importing X.
- Some Mathlib oleans aren't pre-built — match imports to a sibling forge file that compiles.

## BANKING LOCATIONS
- `~/projects/riemann/LEDGER.md` (Mac) — full per-brick detail, §-numbered.
- Olympus `~/rhforge/STATUS.md` — short status line per brick.
- `~/.claude-cobo/projects/-Users-scones/memory/rh-forge-formal-verification.md` — durable memory arc.
- `~/projects/riemann/PROCESS.md` — cumulative decision chain (append, never overwrite).

---

## CURRENT FRONTIER (as of 2026-08-20, §1274 — update as it moves)
Strongest reduction `riemannHypothesis_of_box_weak_limit'` (RHFromBoxLimitClean): confinement + positivity
DISCHARGED; owes only **himz** + **heig**.
- **himz** (no real ζ-zeros in (0,1)) → reduced (§1272/1273/1274) to the single **eta identity**
  `ζ(σ)(1-2^{1-σ})=η(σ)` on (0,1). Positivity η>0 proven (§1273). Remaining: even/odd tsum split + continuation.
- **heig** = the spectral/Weil-Hodge construction = the deep irreducible crux.
Also: analytic sharp zero-free region owes only the log-growth bound (main term §1264 + tail §1268 + assembly
§1269 done; owes the truncated-representation derivation). One-sided RH criterion §1270/1271 unifies lanes:
RH ⟺ keep zeros out of the right half-strip.
