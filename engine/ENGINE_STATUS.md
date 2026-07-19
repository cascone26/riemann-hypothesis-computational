# RH ENGINE — status handoff (2026-07-19)

## STATE: alive on Olympus, and NOW ACTUALLY WORKING (was silently inert)

### The bug that mattered (found + fixed 2026-07-19)
The daemon had run **956 cycles doing nothing**. Two root causes, both fixed:
1. **Generation never parsed.** `call_free` used `ollama run` (the interactive TUI), which
   wrap-redraws lines and duplicates fragments (`sc`→`scipy`, `analo`→`analogous`) — `json.loads`
   choked every time, so the registry never grew past its 4 seeds. **Fix:** `call_free` now hits the
   Ollama **HTTP API** (`localhost:11434/api/generate`) → clean JSON, no TUI artifacts. Verified live.
2. **Seed method code wasn't on Olympus.** The 4 seeds pointed at Mac-only operators → every cycle
   found nothing runnable → idled. **Fix:** synced the seeds + shared dep `critical_eigvec_id.py`.

### New: smoke + free-repair gate on generation
Free qwen writes plausible-but-often-buggy numpy. Now each generated method is **smoke-run (60s)
before it registers**; a crash gets **one free repair pass** (hand the traceback back to the model);
still-broken methods are skipped, not queued. Keeps the loop from clogging with dead code. Auto
methods run under system `python3` (self-contained).

### Confirmed working end-to-end (2026-07-19)
- Generation via FREE Ollama HTTP adds real methods (function_field/selberg/lee_yang/jensen/
  cohomology/spectral…), smoke-gated. Zero Max quota, zero cash.
- Loop runs methods, logs to `outcomes.jsonl`, classifies, escalates on stagnation, refills queue.
- Honest yield: free-generated probes are mostly NULL/ERROR (shallow). The value is **breadth at zero
  cost that never halts**; a real LEAD escalates, and the verify gate blocks false proofs.

### Unchanged safety spine (all intact)
- **API-key HARD GUARD** — `ANTHROPIC_API_KEY/AUTH_TOKEN` scrubbed everywhere; no path bills the key.
- **Adaptive budget** (budget.py) — 0–20 Claude calls/day by week headroom; FREE-only >60%; breaker.
- **Verify gate** (verify_gate.py) — PROOF claim is REPORT until independent Arb/Lean re-check = HANDLE.
- **ntfy** — breakthroughs only, ALL-CAPS, only on a gate-verified proof.
- **Daemon** — `run_forever.sh` on Olympus (heartbeat, halts only on `results/PROOF_FOUND.marker`).

## CONFIRM / RESTART
Check:   ssh olympus 'echo age=$(( $(date +%s) - $(cat ~/projects/riemann/engine/heartbeat) ))s; pgrep -fc run_forever'
Restart: ssh olympus 'cd ~/projects/riemann/engine && nohup bash run_forever.sh >forever.out 2>&1 &'
Watch:   ssh olympus 'tail -f ~/projects/riemann/engine/outcomes.jsonl'

## Files (engine/)
rh_engine.py (loop), generate_methods.py (FREE HTTP generator + smoke/repair gate), budget.py,
verify_gate.py, run_forever.sh (daemon), ENGINE.md (design), methods.json (registry),
outcomes.jsonl (learning log), heartbeat (liveness).
