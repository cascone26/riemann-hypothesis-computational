#!/usr/bin/env python3
"""
t12_roundtable_gate.py -- the non-peer conductor gate for Roundtable-sourced claims.

Why this exists (2026-07-28 design session): multi-agent-debate literature shows symmetric
peer deliberation (N equally-capable seats, e.g. this project's 6-Fable Roundtable) measures up
to 99.2% escalation-of-commitment in some studies, vs ~46% for asymmetric/hierarchical
structures (arXiv:2508.01545). An in-room "skeptic" persona erodes under conformity pressure
over a long session -- majority-conformity effects show up even when the dissenter is right
(arXiv:2509.05396) -- and debate framing can itself destabilize personas rather than stabilize
dissent. The fix isolated-verifier research supports (arXiv:2605.00914, isolated self-correction
beating unguided multi-agent debate) is a checker structurally OUTSIDE the peer group, not just
another seat at the table.

This tool is that checker for this project. It is NOT an LLM persona -- it is a deterministic
gate any Roundtable "this looks novel" / "this is the keystone" claim must pass through, reusing
the SAME honesty spine (core.honesty) and prior-art DB (t01_route_prefilter) every other arsenal
tool already answers to, before the claim is allowed to be banked as more than CLAIMED.

Usage:
    from t12_roundtable_gate import gate_claim
    result = gate_claim({... same 6 boolean constraints as t01, plus lit_check_logged/note ...})
    python3 t12_roundtable_gate.py     # demo on preset roundtable-style claims + selftest
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core.honesty import Result, Tier, guard_text
from core import provenance
from t01_route_prefilter import assess


def gate_claim(claim: dict) -> dict:
    """
    claim: {
      name, description,
      is_global, is_class_generic, is_nonstatistical, is_lower_bound, is_limit_sharp, has_sign_bridge,
      lit_check_logged: bool  -- did the session actually run a literature check for this claim
                                 (t01 prior-art DB and/or engine/lit_recheck.py's LIT_WATCH.md)?
      lit_check_note: str     -- what was checked / against what sources, for the audit trail
    }
    Returns {verdict (from t01), tier (CLAIMED/REPORT), reason, result (honesty.Result)}.
    A claim reaches REPORT (gate-passed, still not proof) ONLY if BOTH:
      (a) t01.assess() verdict == PASSES-FILTER (no DEAD/WALL/KNOWN hit, no failed/unknown constraint)
      (b) lit_check_logged is True with a non-empty note.
    Anything else is held at CLAIMED -- asserted, ungated, not yet earning REPORT tier. This is
    the mandatory checkpoint: peer agreement in the roundtable is not itself evidence.
    """
    v = assess(claim)
    lit_ok = bool(claim.get("lit_check_logged")) and bool((claim.get("lit_check_note") or "").strip())
    passed = (v["verdict"] == "PASSES-FILTER") and lit_ok

    name = claim.get("name", "unnamed")
    if passed:
        tier = Tier.REPORT
        reason = "gate PASSED: no prior-art hit, all six constraints satisfied, lit-check logged."
    else:
        tier = Tier.CLAIMED
        reasons = []
        if v["verdict"] != "PASSES-FILTER":
            reasons.append(f"t01 verdict={v['verdict']}")
        if not lit_ok:
            reasons.append("no logged lit-check (lit_check_logged/lit_check_note missing)")
        reason = "gate HELD at CLAIMED: " + "; ".join(reasons)

    claim_text = guard_text(f"Roundtable claim '{name}': {reason}")
    result = Result(claim=claim_text, tier=tier, method="t12_roundtable_gate + t01 prior-art",
                     tool="t12_roundtable_gate",
                     data={"t01": v, "lit_check_logged": bool(claim.get("lit_check_logged", False))})
    return {"verdict": v["verdict"], "tier": str(tier), "reason": reason, "result": result}


# ---- preset demo claims (modeled on real Roundtable-style claims from LEDGER/PROCESS) ----
PRESETS = [
    {"name": "GORZ-style Jensen ladder rediscovery",
     "description": "closed degrees 3-25 hyperbolicity ladder, claimed novel",
     "is_global": True, "is_class_generic": True, "is_nonstatistical": True,
     "is_lower_bound": True, "is_limit_sharp": False, "has_sign_bridge": False,
     "lit_check_logged": False, "lit_check_note": ""},
    {"name": "Hodge-index arithmetic-surface keystone",
     "description": "extend proven arithmetic Hodge index (Faltings-Hriljac/Yuan-Zhang) past its "
                     "fibration hypothesis to the Spec Z self-square",
     "is_global": True, "is_class_generic": True, "is_nonstatistical": True,
     "is_lower_bound": True, "is_limit_sharp": True, "has_sign_bridge": True,
     "lit_check_logged": True,
     "lit_check_note": "checked Faltings-Hriljac, Yuan-Zhang adelic line bundles, Buium delta_p, "
                        "Burnol Sonine=de Branges, Weil-1964 metaplectic -- see LEDGER 329-332"},
    {"name": "Beurling majorant retry",
     "description": "pointwise majorant dominating the prime comb",
     "is_global": False, "is_class_generic": True, "is_nonstatistical": True,
     "is_lower_bound": True, "is_limit_sharp": False, "has_sign_bridge": False,
     "lit_check_logged": True, "lit_check_note": "checked OBSTRUCTION_MAP r2"},
]


def selftest(verbose=True) -> bool:
    ok = True
    r1 = gate_claim(PRESETS[0])
    if r1["tier"] != "CLAIMED":
        ok = False; print("FAIL: unlit-checked rediscovery should be held at CLAIMED")
    r2 = gate_claim(PRESETS[1])
    if r2["tier"] != "REPORT":
        ok = False; print("FAIL: clean, lit-checked, constraint-passing claim should reach REPORT")
    r3 = gate_claim(PRESETS[2])
    if r3["tier"] != "CLAIMED":
        ok = False; print("FAIL: known-wall route should be held at CLAIMED even if lit-checked")
    for label, r in [("rediscovery", r1), ("keystone", r2), ("beurling", r3)]:
        if verbose:
            print(f"  [{r['tier']:8s}] {label}: {r['reason']}")
    print("t12.selftest:", "PASS" if ok else "FAIL")
    return ok


def main():
    print("=== t12 roundtable gate -- preset claims ===")
    results = []
    for c in PRESETS:
        r = gate_claim(c)
        print(f"\n{c['name']}\n  tier: {r['tier']}  verdict: {r['verdict']}\n  reason: {r['reason']}")
        results.append(r["result"])
    print()
    ok = selftest(verbose=False)
    path = provenance.stamp("t12_roundtable_gate", results, inputs={"n_presets": len(PRESETS)},
                            extra={"selftest_pass": ok})
    print("run record:", path)
    return ok


if __name__ == "__main__":
    main()
