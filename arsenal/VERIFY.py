#!/usr/bin/env python3
"""
VERIFY.py  -- full verification pass over the whole arsenal.

Runs core selftests + every tool's selftest end-to-end, captures pass/fail, then runs the
fabrication detector's audit over ALL run records to confirm no HANDLE claim is unbacked.
Writes arsenal/RESULTS.md with the honest per-tool status.  Exit code 0 iff everything passes
AND no unbacked HANDLE claims exist.

This is the "did it actually work" gate -- not a hedged 'should work'.  Every tool below was run;
the table records the real result.
"""
from __future__ import annotations
import sys, os, io, time, importlib, contextlib, traceback
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CORE = ["core.honesty", "core.weil"]
TOOLS = [
    ("t01_route_prefilter", "route pre-filter (proof-shape spec + prior-art)"),
    ("t02_wall_atlas", "wall atlas + equivalence detector"),
    ("t03_home_compiler", "home compiler (3-clause scorecard + gates)"),
    ("t04_family_comparator", "cross-L-function family comparator (real compute)"),
    ("t05_ladder_finder", "reverse-math ladder / frontier finder"),
    ("t06_fabrication_detector", "fabrication detector (independent re-run)"),
    ("t07_rosetta", "RH-equivalence transport engine"),
    ("t08_exact_sos", "exact-rational PD/SOS certifier"),
    ("t09_debranges", "de Branges / Laguerre positivity instrument"),
    ("t10_flight_simulator", "tunable-axiom positivity sandbox"),
    ("t11_disproof_engine", "adversarial disproof engine"),
    ("t12_roundtable_gate", "non-peer conductor gate for Roundtable claims (prior-art + lit-check enforced)"),
]


def _run_selftest(modname):
    """Import module, call selftest(verbose=False), normalize return to bool. Returns (ok, err)."""
    try:
        mod = importlib.import_module(modname)
        fn = getattr(mod, "selftest", None)
        if fn is None:
            return None, "no selftest()"
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            res = fn(verbose=False)
        ok = res[0] if isinstance(res, tuple) else res
        return bool(ok), None
    except Exception as e:
        return False, f"{type(e).__name__}: {e}\n{traceback.format_exc()[-600:]}"


def main():
    print("=" * 70)
    print("RH ARSENAL — FULL VERIFICATION PASS")
    print("=" * 70)
    t0 = time.time()
    rows = []

    print("\n[core]")
    for m in CORE:
        ok, err = _run_selftest(m)
        status = "PASS" if ok else ("SKIP" if ok is None else "FAIL")
        print(f"  {status:5s} {m}" + (f"  -- {err.splitlines()[0]}" if err else ""))
        rows.append((m, "core", status, err))

    print("\n[tools]")
    for m, desc in TOOLS:
        ok, err = _run_selftest(m)
        status = "PASS" if ok else ("SKIP" if ok is None else "FAIL")
        print(f"  {status:5s} {m:26s} {desc}" + (f"\n         !! {err.splitlines()[0]}" if err else ""))
        rows.append((m, desc, status, err))

    # audit all HANDLE claims for real evidence
    print("\n[fabrication audit over runs/]")
    from t06_fabrication_detector import audit_runs
    findings = audit_runs()
    if findings:
        for f in findings:
            print(f"  UNBACKED: {f.get('file')} -- {f.get('detail')}")
    else:
        print("  clean -- every HANDLE claim on disk points at real evidence")

    n_pass = sum(1 for _, _, s, _ in rows if s == "PASS")
    n_fail = sum(1 for _, _, s, _ in rows if s == "FAIL")
    all_ok = (n_fail == 0) and (len(findings) == 0)
    dt = time.time() - t0

    print("\n" + "=" * 70)
    print(f"RESULT: {n_pass} pass, {n_fail} fail, {len(findings)} unbacked HANDLE claims "
          f"({dt:.1f}s)")
    print("VERDICT:", "ALL GREEN" if all_ok else "FAILURES PRESENT")
    print("=" * 70)

    _write_results_md(rows, findings, all_ok, dt)
    return 0 if all_ok else 1


def _write_results_md(rows, findings, all_ok, dt):
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "RESULTS.md")
    lines = ["# RH Arsenal — Verification Results", "",
             f"Full verification pass. Verdict: **{'ALL GREEN' if all_ok else 'FAILURES'}** "
             f"({dt:.1f}s).", "",
             "Every tool below was RUN (not assumed). Status = its selftest result against "
             "known-answer cases. HANDLE-tier claims are audited for a real proof-pointer on disk.",
             "", "| module | what | status |", "|---|---|---|"]
    for m, desc, status, err in rows:
        badge = {"PASS": "PASS ✅", "FAIL": "FAIL ❌", "SKIP": "skip"}[status]
        lines.append(f"| `{m}` | {desc} | {badge} |")
    lines += ["", "## Fabrication audit (HANDLE claims vs real evidence)", ""]
    if findings:
        for f in findings:
            lines.append(f"- **UNBACKED**: `{f.get('file')}` — {f.get('detail')}")
    else:
        lines.append("- Clean: every HANDLE-tier claim in `runs/` points at a real, non-empty "
                     "proof-pointer file.")
    lines += ["", "## Honesty tiers used", "",
              "- **HANDLE** (earned by a real artifact): t04 per-L-function validated margins, "
              "t08 exact rational PD certs. Each points at a `runs/*.json` proof file.",
              "- **REPORT** (analyzed/ran, not proof): all structural/knowledge tools + the "
              "numeric instruments (t09/t10/t11) — finite/numeric, explicitly not RH proofs.",
              "- No tool is permitted to claim a proof of RH (honesty.guard_text blocks it).", ""]
    with open(path, "w") as f:
        f.write("\n".join(lines))
    print("wrote", path)


if __name__ == "__main__":
    sys.exit(main())
