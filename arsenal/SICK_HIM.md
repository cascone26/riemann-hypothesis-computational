# SICK_HIM — adversarial verification protocol for the RH arsenal

"SICK HIM" = the Foundry's honesty gate, applied to this toolkit: **assume every claim is a lie
until an INDEPENDENT check says otherwise.** A green selftest is not enough — a selftest can be a
rubber stamp. This pass attacks the arsenal to catch false certifications, overclaims, and math
errors before any of it is trusted.

Motivation: the project has twice been burned by fabricated "CERTIFIED" claims that only fell to
independent re-run. So nothing here is believed on the strength of its own report.

## The seven adversarial checks (`SICK_HIM.py`)

1. **Overclaim sweep** — grep all tool source AND all `runs/*.json` emitted claims for banned
   RH-proof phrasings ("proven RH", "certified RH", …). Any hit = fail.
2. **HANDLE evidence deep-check** — for every HANDLE-tier claim on disk, open its proof-pointer
   file and confirm it exists, is non-empty, and actually contains the asserted content (not a
   placeholder). A HANDLE with hollow evidence = fail.
3. **Triple-independent Weil check** — recompute the Weil margin a THIRD way (raised precision,
   shifted quadrature grid, independent zero source) and require agreement with both the numpy
   engine and the mpmath path. Two agreeing paths could share a bug; a third breaks ties.
4. **Exact-cert adversarial (t08)** — feed near-singular matrices (det = ±1e-18), a PSD-not-PD
   matrix (det = 0), and a non-symmetric matrix; require the exact-rational verdict to be correct
   in every case (cross-checked against the true rational det). Proves float can't fool it.
5. **Mutation probe (the key one)** — inject a SIGN FLIP into the Weil engine at runtime and
   confirm `weil.selftest` THEN FAILS. If the selftest still passes on a broken engine, the
   selftest is a rubber stamp and everything downstream is suspect. This proves the tests bite.
6. **Guard fuzz** — throw a dozen overclaim phrasings at `honesty.guard_text` and confirm each is
   blocked; confirm a HANDLE with a missing evidence file is refused at construction.
7. **Independent-code audit** — confirm `t06.independent_weil_margin` shares NO computation path
   with `core.weil` (both must be able to disagree), so their agreement is meaningful.

## Verdict rule
Clean iff: 0 overclaims, 0 hollow HANDLE claims, triple-check agrees, all t08 adversarial verdicts
correct, the mutation probe DOES break the selftest, and all guard fuzz blocked. Anything else is
reported loudly and the arsenal is NOT trusted until fixed.

Plus: three independent skeptic agents (math errors / overclaims / fraud-resistance) review in
parallel; the control session re-verifies each finding rather than trusting the reports.
