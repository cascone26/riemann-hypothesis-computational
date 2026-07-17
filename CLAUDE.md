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
