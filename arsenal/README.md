# RH Arsenal

Eleven RH-analysis tools built to be **different in kind** from "more numerics / more zeros" — each
exploits an asymmetric asset (cross-domain synthesis, tireless bookkeeping, adversarial self-check)
rather than competing on raw analytic skill. Built 2026-07-26. Every tool is honestly tiered and
survives an adversarial anti-fabrication pass.

## Run it

```bash
cd ~/projects/riemann/arsenal
python3 VERIFY.py       # runs all 13 selftests + audits every HANDLE claim  -> RESULTS.md
python3 SICK_HIM.py     # adversarial pass: overclaims, mutation probes, fraud attempts
python3 t04_family_comparator.py   # (or any t01..t11) to run one tool
```

## The honesty layer (why this can't quietly lie)

`core/honesty.py` enforces three tiers, never conflated:
- **CLAIMED** — asserted, no evidence. **REPORT** — ran/analyzed, not proof. **HANDLE** — earned by
  a real artifact on disk. `assert_handle` refuses HANDLE if the proof-pointer file is absent; a
  banned-phrase guard blocks any "proven/certified RH" overclaim at construction time.
- Only **t04** (validated L-margins) and **t08** (exact rational PD certs) earn HANDLE; everything
  else is REPORT (finite/numeric, explicitly not an RH proof).

## The tools

| # | tool | tier | what it does |
|---|------|------|--------------|
| 1 | `t01_route_prefilter` | REPORT | scores a proposed approach vs the proof-shape spec (§4 i–vi) + prior-art DB; flags dead/known routes before a session is spent |
| 2 | `t02_wall_atlas` | REPORT | graph of routes→obstructions→reductions; equivalence detector shows 8 classical routes collapse to one Wall Invariant |
| 3 | `t03_home_compiler` | REPORT | scores candidate geometric "homes" vs the 3-clause scorecard + 2 gates; pins the missing signature/pairing keystone |
| 4 | `t04_family_comparator` | **HANDLE** | real Weil margin across ζ + 6 Dirichlet L-functions, each validated vs its own zeros; positivity is family-universal |
| 5 | `t05_ladder_finder` | REPORT | reverse-math lattice of proven→RH; names the nearest open stepping-stone rung (support_C = the continuum target) |
| 6 | `t06_fabrication_detector` | REPORT/HANDLE | the anti-fraud tool: re-derives claims by an independent path, audits HANDLE evidence, deep-reverifies stored certs |
| 7 | `t07_rosetta` | REPORT | RH-equivalence transport engine; moves a partial result between the ~10 equivalent forms |
| 8 | `t08_exact_sos` | **HANDLE** | exact-rational LDL^T PD certs over ℚ (no floats); reproduces the moment route's death exactly |
| 9 | `t09_debranges` | REPORT | de Branges/Hermite-Biehler positivity via the Laguerre inequality on ξ; ties de Branges to the Jensen ladder |
| 10 | `t10_flight_simulator` | REPORT | tunable-axiom sandbox; maps the forcing boundary and shows positivity tightening as t→limit |
| 11 | `t11_disproof_engine` | REPORT | adversary: bends prime axioms, reports which bend breaks Weil positivity (PNT/Euler are the load-bearing ones) |

## Shared engine

`core/weil.py` — the Riemann-Weil explicit-formula functional A+B−P, self-validating: it computes
both the analytic side and the zero side and reports the closing residual (~1e-12 for ζ). That
residual is the engine's own correctness certificate. Parametric, so t04/t10/t11 can bend the axioms.

## Anti-fabrication (SICK_HIM)

`SICK_HIM.py` runs 8 adversarial checks: overclaim sweep, HANDLE evidence deep-check, a **triple-
independent** Weil recompute, exact-cert adversarial matrices, a **mutation probe** (injects a sign
flip and confirms the selftest then FAILS — proving the tests aren't rubber stamps), guard fuzz,
code-independence audit, and audit-hardening (dir/symlink/junk evidence + a fabricated cert caught
by re-derivation). Three independent skeptic agents also reviewed math / overclaims / fraud-
resistance; their findings (a dead-but-wrong helper, an audit bypass) were verified and fixed.

## Honest scope

Nothing here proves RH. These are instruments and gates that (a) stop wasted sessions, (b) attack
the confirmed keystone (the missing definite-signature pairing), and (c) make our own claims
impossible to fake without an independent check catching it. See `BUILD_PLAN.md`, `RESULTS.md`,
`SICK_HIM.md`.
