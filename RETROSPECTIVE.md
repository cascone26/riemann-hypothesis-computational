# RH Forge — Process Retrospective & Speedup Roadmap
_Written 2026-08-21, at Jacob's request, ahead of riemann9. Based on direct data from the §1352→§1630 run (279 EXIT0 bricks this cycle)._

This is an honest audit of **what actually slows us down** and **the speedups we haven't been using** — including
things that are genuinely available right now that we just haven't reached for. Ranked by leverage.

---

## Part 1 — What's actually slowing us down (measured, this session)

### 1. Mathlib API name-churn is the #1 time sink (~30–40% of iterations)
Nearly every brick cost 1–3 recompile round-trips **not on the math, but on lemma names/signatures that moved or
that I guessed wrong.** Concrete failures this session alone:
`le_or_lt` (gone → `le_total`/`by_cases`, hit 3×), `tsum_le_tsum` (→ `Summable.tsum_le_tsum`),
`Function.update_same` (→ `update_self`), `div_le_div_of_nonneg_right` (wrong args),
`Ioo_mem_nhdsWithin_Ioi` (→ `Ioo_mem_nhdsGT`), `norm_term_eq` (has a hidden `if n=0` guard),
`ContinuousAt.comp` / `.div` (beta-redex unification needs an explicit target-type `have`),
`const_cpow` (not reachable by `fun_prop`). Each = scp + compile (~2–4s) + my reasoning to re-guess.
**This is pure, avoidable drag.** The math was usually right the first time; the *names* weren't.

### 2. Per-brick banking overhead (~7 manual steps each)
Every single brick: write file → scp to Olympus → `lean` compile → check `#print axioms` → `lean -o` olean →
`scp` mirror back to Mac → append LEDGER → append STATUS → append memory. ~7 steps, ~half a dozen tool calls,
for a 1-theorem file. Over 279 bricks that's ~1,900 mechanical steps.

### 3. Everything built serially, by hand
One brick at a time. Many bricks are **independent** (different files, no shared dependency) yet I compile them
one after another. And — the big one — **I hand-authored all 279**, while a working AI prover loop sat idle (see
Part 2.1).

### 4. Silent bit-rot (found today, would have stayed hidden)
Today's verification sweep found **4 files that no longer compile** (Mathlib renamed lemmas under them, e.g.
`Real.summable_one_div_nat_rpow`) and **13 "done" modules with no persisted `.olean`** — including the zero-free
capstones §1575/§1577. All dead-leaf or now-restored, tower intact — but we only know that because we *looked*.
Without a standing CI, "done" silently decays.

### 5. Context accumulation (single-session constraint)
`rh-truncated-rep-chain.md` is 243 lines; the LEDGER is enormous. Recall and re-orientation cost grows with the
run. (Single-session is Jacob's confirmed preference — so the fix is leaner memory, not multi-session.)

---

## Part 2 — Speedups we're NOT using but genuinely could (the "you don't think of" list)

### 1. ★ The AI prover loop is ALIVE and IDLE — biggest untapped lever
A Qwen3-32B-AWQ vLLM server is **running on Olympus's 3090s right now** (`localhost:8000`, 0% GPU util =
idle), with `pipeline/prove_loop2_v{3,4,6}.py` generate→verify loops already built. **I've been hand-writing
every brick while an automated prover sits unused.** Even if it only lands the easy 40% (routine
integrability/continuity/algebra bricks), that's a throughput multiplier — and our new **axiom gate verifies its
output** so it can't fool us. Caveat: the server's `max_model_len` is only 2048 tokens (tiny) — fine for small
bricks, needs a bump for big ones. **Action: wire prove_loop → axiom_gate as a first-draft generator; I curate.**

### 2. ★ LeanSearchClient / `exact?` / `apply?` are installed but unused
Mathlib's `LeanSearchClient` package **is present** in our toolchain. `exact?`, `apply?`, and type-directed
search kill the #1 bottleneck (name-churn) directly — query by goal type instead of guess→compile→fix.
**Action: run `exact?`/`apply?` in-place before hand-writing a term, and keep a verified-name cheatsheet.**

### 3. A Lean **LSP server** for instant feedback instead of scp+compile round-trips
Right now every check is a full `lean` invocation over SSH (~2–4s + my parse of the error). A persistent
`lean --server` (or editing against the cached Mathlib) gives real-time hover/goal-state/error — collapsing the
iteration loop from seconds-per-attempt to instant. **This is the single biggest iteration-latency win.**

### 4. `bank.sh` — one command for the whole 7-step ritual
A single `forge-brick <file> <Module> "<ledger line>"` that scp's, compiles, axiom-audits, mirrors, and appends
LEDGER/STATUS. Turns ~7 steps into 1. (Building this today.)

### 5. Batch/parallel compilation of independent bricks
Draft N independent bricks, `scp` together, compile concurrently (`&`/parallel). Wall-clock ≈ slowest single
brick instead of the sum.

### 6. A dependency-DAG + true-frontier report
We track progress by hand-numbered §. An auto-computed import graph shows the real *leaf* theorems (the actual
frontier), flags orphans/dead-leaves (would have caught today's 9 dead files instantly), and makes circular
reasoning impossible to hide. **Also tells us the highest-leverage next brick instead of a random tail estimate.**

---

## Part 3 — The honest ceiling (what tooling can and cannot do)

None of the above proves RH. The actual crux — extending the arithmetic Hodge-index positivity to
Spec ℤ ×_F₁ Spec ℤ (Connes–Consani), or a non-circular spectral realization (heig / Connes) — is **research-scale
new mathematics the forge cannot invent**, and the forge memory names *fabrication* (frontier math that looks
right and isn't) as the true blocker. So tooling's job is precisely scoped:
- **maximize verified-scaffolding throughput** (prover loop, LSP, bank.sh, batching),
- **make it impossible to fool ourselves** (axiom gate ✓, novelty/Report-Handle gate, canary, CI),
- **always work the highest-leverage next thing** (frontier DAG).

That is the real answer to "make our work always the best": not a cleverer tail-estimate, but a machine that
generates fast, verifies un-foolably, and never lets a broken or fabricated brick pass as done.

---

## Ranked build order (recommendation)
1. ✅ **Axiom gate** — built today; already caught 9 hidden holes + 4 bit-rotted files.
2. **`exact?`/LeanSearchClient into the workflow** + a verified-name cheatsheet — kills the #1 bottleneck, zero build cost.
3. **Reactivate the vLLM prover loop → axiom gate** — the throughput multiplier that's already 90% built and idle.
4. **`bank.sh`** — kills the #2 bottleneck (mechanical overhead). (Building today.)
5. **Standing CI** (`axiom_gate.sh` on a cron) — prevents silent bit-rot recurring.
6. **Novelty / Report-Handle gate** — the direct anti-fabrication lever (the failure mode that most threatens quality).
7. **Dependency-DAG + frontier report** — leverage-targeting + orphan detection.
8. **Lean LSP** — the deepest iteration-latency win (higher setup cost).
