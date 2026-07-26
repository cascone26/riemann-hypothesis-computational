# RH Arsenal — Build Plan (2026-07-26, overnight autonomous build)

Goal: build the Tier S + A + B tools from the 20-tool prism ranking. Every tool honestly
scoped, every claim carrying a proof-pointer, and a full verification + adversarial SICK_HIM
pass so nothing gets a false certification.

## Non-negotiable discipline (why this arsenal exists at all)
The project's two worst repeated failures are (1) rediscovering published work and (2)
fabricating our own "CERTIFIED" claims (2 of 3 parallel tracks faked certs, caught only by
independent re-run). So the arsenal is built AROUND an honesty layer:

- **Three tiers, never conflated** (Charisma Report/Handle + our proof-pointer rule):
  - `CLAIMED`  — asserted, no evidence yet. Default for anything unrun.
  - `REPORT`   — describable/analyzable now (a scorer ran, a structural check passed). Real work, not proof.
  - `HANDLE`   — earned ONLY by a real run producing a real artifact (a rigorous Arb cert, an
    exact rational LDL^T with all pivots > 0, a re-derivation that agrees). A claim cannot be
    stamped HANDLE without a proof-pointer file that exists on disk.
- Every result is a `Result` object carrying tier + proof_pointer + method. `assert_handle`
  refuses to mark HANDLE if the evidence file is absent.
- No tool prints "proven RH" / "certified RH". Finite/numeric results say so in the same breath.

## The engine anchors (search-before-build — these already exist and are validated)
- Weil form A(g)+B(g)−P(g), Ω(r)=Re ψ(¼+ir/2)−log π, primes Σ Λ(k)k^{−1/2}, boundary W_{0,2}
  — from `operators/connes_decisive_v3.py`, `operators/auto_hodge_index_signature.py`.
- Wall structure + proof-shape spec (i–vi) — `RH_WALL_ANATOMY.md` §4, `OBSTRUCTION_MAP.md`.
- Disproof/Beurling-comb idea already seeded in `operators/beurling_adversary.py` (we generalize).
- Home/Hodge-index scoring seeded in `operators/auto_hodge_index_*.py` (we wrap + score).

## The 11 tools (honest v1 scope each)

### Tier S
1. **t01 route_prefilter** — rule engine over the proof-shape spec (§4 i–vi) + Wall Invariant +
   a known-published-results DB (Dimitrov–Lucas, GORTTW, obstruction routes 1–7). Input: a
   structured description of a proposed approach. Output: which constraints it violates, whether
   it's sign-blind (can't beat the explicit formula), and nearest known prior art. REPORT-tier
   classifier. Value: stops us investing in dead/known routes.
2. **t02 wall_atlas** — parse routes→obstructions→reductions into a graph (networkx). Equivalence
   detector v1 = find routes that reduce to the same obstruction node (the Wall Invariant made
   mechanical). REPORT-tier structural tool over already-verified data.
3. **t03 home_compiler** — scorer, NOT a synthesizer. Encodes the 3-clause Roundtable scorecard
   (prior / generative / signature-bearing) + two hard gates (reproduce explicit formula as a
   Lefschetz trace; carry a definite-signature pairing). Scores candidate "homes" (C-C spectral,
   C-C geometric RR, de Branges, Nyman-Beurling, function-field analogue). REPORT-tier.

### Tier A
4. **t04 family_comparator** — REAL compute: Weil positivity margin across a family (ζ, Dirichlet
   L(χ) for a few χ) via the shared engine; reports what's invariant across the family. Finite-N,
   HANDLE-tier for each computed margin (proof-pointer = run JSON), REPORT-tier for the "invariant".
5. **t05 ladder_finder** — logical lattice of statements between what we've proven (n≤9.8M,
   pair-positivity, function-field RH) and RH, ordered by strength; flags candidate stepping-stones
   with a not-yet-crossed gap. REPORT-tier reasoning scaffold.
6. **t06 fabrication_detector** — the anti-fraud tool. Takes a claim + its proof-pointer (a result
   JSON / a script), RE-RUNS it independently and/or cross-checks numeric values against an
   independent method, flags disagreement or missing evidence. HANDLE-tier verdicts (it actually
   re-executes). Directly targets failure mode (2).

### Tier B
7. **t07 rosetta** — the ~20 RH-equivalents as a graph with transport annotations; given a partial
   result in one form, report its image/strength in the others. REPORT-tier knowledge engine.
8. **t08 exact_sos** — sympy exact-rational LDL^T / SOS attempt on small tau matrices (the certified
   first-prime 4×4 window). If all pivots > 0 over ℚ, that's a genuine EXACT PD certificate for that
   finite matrix (HANDLE-tier, proof-pointer = the rational pivots). Honest: finite matrix, not RH.
9. **t09 debranges** — build the canonical-system Hamiltonian H(t) from ξ and test the structure
   function |E| condition numerically at small scale. REPORT-tier numeric instrument.
10. **t10 flight_simulator** — tunable-axiom sandbox on the shared engine: dial arch factor / prime
    distribution / FE symmetry, watch the positivity margin, map the forcing boundary. REPORT-tier.
11. **t11 disproof_engine** — generalizes beurling_adversary: search bent-axiom L-function families
    for an off-line-zero-carrying object, report which axiom is forced to break. REPORT-tier.

## Deliverables
- `arsenal/core/{honesty,provenance,weil}.py` — trust foundation + shared engine.
- `arsenal/t01..t11_*.py` — the 11 tools, each with a `__main__` demo + a `selftest()`.
- `arsenal/VERIFY.py` — runs all, captures to runs/, asserts proof-pointers exist.
- `arsenal/SICK_HIM.py` + `SICK_HIM.md` — adversarial re-verification harness + protocol.
- `arsenal/RESULTS.md` — honest per-tool status (tier + proof-pointer) after verification.
