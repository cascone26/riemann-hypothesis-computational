# RH Arsenal — Verification Results

Full verification pass. Verdict: **ALL GREEN** (44.4s).

Every tool below was RUN (not assumed). Status = its selftest result against known-answer cases. HANDLE-tier claims are audited for a real proof-pointer on disk.

| module | what | status |
|---|---|---|
| `core.honesty` | core | PASS ✅ |
| `core.weil` | core | PASS ✅ |
| `t01_route_prefilter` | route pre-filter (proof-shape spec + prior-art) | PASS ✅ |
| `t02_wall_atlas` | wall atlas + equivalence detector | PASS ✅ |
| `t03_home_compiler` | home compiler (3-clause scorecard + gates) | PASS ✅ |
| `t04_family_comparator` | cross-L-function family comparator (real compute) | PASS ✅ |
| `t05_ladder_finder` | reverse-math ladder / frontier finder | PASS ✅ |
| `t06_fabrication_detector` | fabrication detector (independent re-run) | PASS ✅ |
| `t07_rosetta` | RH-equivalence transport engine | PASS ✅ |
| `t08_exact_sos` | exact-rational PD/SOS certifier | PASS ✅ |
| `t09_debranges` | de Branges / Laguerre positivity instrument | PASS ✅ |
| `t10_flight_simulator` | tunable-axiom positivity sandbox | PASS ✅ |
| `t11_disproof_engine` | adversarial disproof engine | PASS ✅ |
| `t12_roundtable_gate` | non-peer conductor gate for Roundtable claims (prior-art + lit-check enforced) | PASS ✅ |

## Fabrication audit (HANDLE claims vs real evidence)

- Clean: every HANDLE-tier claim in `runs/` points at a real, non-empty proof-pointer file.

## Honesty tiers used

- **HANDLE** (earned by a real artifact): t04 per-L-function validated margins, t08 exact rational PD certs. Each points at a `runs/*.json` proof file.
- **REPORT** (analyzed/ran, not proof): all structural/knowledge tools + the numeric instruments (t09/t10/t11) — finite/numeric, explicitly not RH proofs.
- No tool is permitted to claim a proof of RH (honesty.guard_text blocks it).
