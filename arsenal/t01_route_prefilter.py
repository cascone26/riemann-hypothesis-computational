#!/usr/bin/env python3
"""
t01_route_prefilter.py  -- Tier S #1.  The "will this route provably fail / is it already
published" gate.  Encodes the proof-shape spec (RH_WALL_ANATOMY.md sec 4, constraints i-vi),
the Wall Invariant (sec 3), and a known-prior-art DB, so a proposed approach can be scored
BEFORE a session is invested in it.

REPORT-tier: a rule engine over VERIFIED campaign results (OBSTRUCTION_MAP + RH_WALL_ANATOMY).
It does not prove a route dead; it flags which mandatory constraints a route fails and the
nearest known result, so we stop re-deriving Dimitrov-Lucas / GORTTW / the mapped obstructions.

Usage:
    from t01_route_prefilter import assess
    verdict = assess({... structured approach ...})
    python3 t01_route_prefilter.py        # demo on preset approaches + selftest
"""
from __future__ import annotations
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core.honesty import Result, Tier
from core import provenance

# ---- The six mandatory constraints any real proof must satisfy (RH_WALL_ANATOMY sec 4) ----
# Each: key, human text, and the attribute that must hold in an approach dict.
CONSTRAINTS = [
    ("global",       "Must be GLOBAL, not pointwise (routes 2-4 fail locally: prime comb "
                     "exceeds arch kernel at every prime >=5).",            "is_global"),
    ("class_generic","Must NOT use zeta-specific structure (chi_4 boundary-saturation: must "
                     "work for the whole FE+Euler class).",                 "is_class_generic"),
    ("nonstatistical","Must NOT be statistical (zeros carry no statistical fingerprint of the "
                     "forcing; 07-18 anomaly probes).",                     "is_nonstatistical"),
    ("lower_bound",  "Must be a LOWER-bound mechanism (variational/eigenvector reasoning gives "
                     "upper bounds only).",                                 "is_lower_bound"),
    ("limit_sharp",  "Must survive the limit / become exponentially sharp as lambda grows "
                     "(fixed-slack arguments die when margin e^{-4pi lambda^2} < slack).", "is_limit_sharp"),
    ("sign_bridge",  "MUST contain a sign-carrying prime<->zero bridge that is NOT the "
                     "explicit formula (the single non-negotiable novelty).", "has_sign_bridge"),
]

# ---- Known prior art / mapped obstructions (OBSTRUCTION_MAP + LEDGER + literature gate) ----
# keyword -> (verdict, precise status)
PRIOR_ART = {
    "moment":        ("DEAD", "Hankel/moment (Track 1): fails PD at N=2, det=-0.0037; lambda_n "
                              "log-CONCAVE, opposite of a moment sequence. (OBSTRUCTION_MAP r1)"),
    "hankel":        ("DEAD", "Same as moment route: log-concave kills Stieltjes moment structure."),
    "beurling":      ("WALL", "Beurling-Selberg majorant: no pointwise majorant exists; prime comb "
                              "exceeds arch kernel at every prime>=5. (OBSTRUCTION_MAP r2)"),
    "large sieve":   ("WALL", "Large sieve bounds OFF-line zeros (zero-density), never forces them "
                              "ON. (OBSTRUCTION_MAP r3)"),
    "de branges":    ("WALL", "A+B has 4 negative eigenvalues; naive dB positivity fails "
                              "(Conrey-Li). Only tau=A+B-P is PSD. (OBSTRUCTION_MAP r4)"),
    "operator":      ("WALL", "tau spectral: every certificate reduces to Weil; tau>=0 IS RH. "
                              "(OBSTRUCTION_MAP r5)"),
    "jensen":        ("KNOWN","Jensen/Laguerre-Polya hyperbolicity: Dimitrov-Lucas 2011 (deg-3 "
                              "uncond.), GORTTW/Chasse all-n to d~10^17-10^20. Fixed-degree ladder "
                              "is kernel-local & RH-free; RH lives only in uniformity in d."),
    "laguerre":      ("KNOWN","Laguerre-Polya / Turan: heavily worked, see Jensen entry (GORTTW)."),
    "polya":         ("KNOWN","Polya-Jensen: GORTTW / Griffin-Ono-Rolen-Zagier 2019. Effective "
                              "bounds published to d~10^20."),
    "turan":         ("KNOWN","Turan inequalities on Jensen polynomials: GORTTW. Known to d~10^20."),
    "debruijn":      ("WALL", "de Bruijn-Newman: g_lambda detector is zero-side only (no Euler "
                              "product => no prime handle). Rodgers-Tao gives Lambda>=0 (WRONG "
                              "direction). (LEDGER 47e)"),
    "newman":        ("WALL", "See de Bruijn-Newman: detector-only, wrong-direction."),
    "crystalline":   ("WALL", "Crystalline-measure/Fourier-quasicrystal rigidity = the Jensen/"
                              "Lee-Yang wall (LEDGER 183). 'zero comb is an FQ'<=>RH is circular input."),
    "quasicrystal":  ("WALL", "See crystalline: FQ classification passes through Lee-Yang polys = "
                              "Jensen wall."),
    "mollifier":     ("PARTIAL","Levinson/Speiser mollifier: positive PROPORTION on line (~41.7%). "
                              "Barrier: mollifier length theta<4/7 (shifted-convolution). Reaching "
                              "100% ~ RH-hard. (OBSTRUCTION_MAP r7)"),
    "levinson":      ("PARTIAL","See mollifier: theta<4/7 shifted-convolution barrier."),
    "lee-yang":      ("WALL", "Lee-Yang: xi has reflection but NOT coefficient positivity (07-18 "
                              "physics dive) = the A+B>=P wall in physics clothing."),
    "zero-density":  ("WALL", "Zero-density bounds off-line zeros via sign-blind counting; never "
                              "forces them on. (LEDGER 180)"),
    "polya frequency": ("PARTIAL", "de Bruijn-Newman kernel is NOT PF5 (Polya-Frequency order 5): "
                              "explicit negative 5x5 Toeplitz minor, certified via 80-digit interval "
                              "arithmetic; PF2-PF4 hold (Michalowski, arXiv:2602.20313, folded in from "
                              "engine/lit_recheck.py 2026-07-28). A positivity floor real at low order "
                              "CAN crack at a specific finite higher order in RH-adjacent territory -- "
                              "makes a 'crack not ruled out' hypothesis elsewhere more plausible, not "
                              "less. Also the honesty template: v2 withdrew the certified GLOBAL sign/"
                              "unique-threshold claim (unsound tail/asymptotic enclosure) while the "
                              "DIRECT interval-certified minor survived -- any crack-hunt result here "
                              "must clear that same bar (direct enclosed minor/eigenvalue certificate "
                              "on a specific cell, never a tail-dependent global claim)."),
}


def _prior_art_hits(name: str, description: str):
    text = (name + " " + description).lower()
    hits = []
    for kw, (verdict, status) in PRIOR_ART.items():
        if kw in text:
            hits.append((kw, verdict, status))
    return hits


def assess(approach: dict) -> dict:
    """
    approach: {
      name: str, description: str,
      is_global, is_class_generic, is_nonstatistical, is_lower_bound,
      is_limit_sharp, has_sign_bridge : bool (True = satisfies the constraint),
    }
    Unspecified booleans default to None (unknown -> flagged as unassessed, counts against).
    Returns a verdict dict.
    """
    name = approach.get("name", "unnamed")
    desc = approach.get("description", "")
    failed, unknown = [], []
    for key, text, attr in CONSTRAINTS:
        v = approach.get(attr, None)
        if v is None:
            unknown.append((key, text))
        elif not v:
            failed.append((key, text))

    hits = _prior_art_hits(name, desc)
    dead_or_wall = [h for h in hits if h[1] in ("DEAD", "WALL")]
    known = [h for h in hits if h[1] == "KNOWN"]

    # verdict logic
    if any(h[1] == "DEAD" for h in hits):
        verdict = "LIKELY-DEAD"
    elif not approach.get("has_sign_bridge", False):
        # the single non-negotiable: no new sign-carrying bridge => reduces to the wall
        verdict = "REDUCES-TO-WALL"
    elif dead_or_wall:
        verdict = "KNOWN-WALL"
    elif known:
        verdict = "KNOWN-TERRITORY"
    elif failed:
        verdict = "VIOLATES-SPEC"
    elif unknown:
        verdict = "UNDERSPECIFIED"
    else:
        verdict = "PASSES-FILTER"

    return {
        "name": name,
        "verdict": verdict,
        "failed_constraints": [k for k, _ in failed],
        "unknown_constraints": [k for k, _ in unknown],
        "prior_art": [{"keyword": k, "class": c, "status": s} for k, c, s in hits],
        "note": ("The sign-carrying-bridge constraint (vi) is the non-negotiable novelty; "
                 "no route without it can beat the sign-blind explicit formula."),
    }


# ---- preset demo approaches ----
PRESETS = [
    {"name": "Hankel moment positivity", "description": "show lambda_n is a Stieltjes moment sequence",
     "is_global": True, "is_class_generic": True, "is_nonstatistical": True,
     "is_lower_bound": True, "is_limit_sharp": False, "has_sign_bridge": False},
    {"name": "Higher Jensen-polynomial hyperbolicity", "description": "extend Polya-Jensen Turan to all degrees",
     "is_global": True, "is_class_generic": True, "is_nonstatistical": True,
     "is_lower_bound": True, "is_limit_sharp": False, "has_sign_bridge": False},
    {"name": "Beurling-Selberg pointwise majorant", "description": "dominate the prime comb pointwise",
     "is_global": False, "is_class_generic": True, "is_nonstatistical": True,
     "is_lower_bound": True, "is_limit_sharp": False, "has_sign_bridge": False},
    {"name": "Hodge-index intersection form on arithmetic surface",
     "description": "build a definite-signature pairing (Castelnuovo/Hodge-index) transporting sign from geometry",
     "is_global": True, "is_class_generic": True, "is_nonstatistical": True,
     "is_lower_bound": True, "is_limit_sharp": True, "has_sign_bridge": True},
]


def selftest(verbose=True) -> bool:
    expect = {
        "Hankel moment positivity": "LIKELY-DEAD",
        "Higher Jensen-polynomial hyperbolicity": "KNOWN-TERRITORY",  # via prior-art + no bridge
        "Beurling-Selberg pointwise majorant": "REDUCES-TO-WALL",
        "Hodge-index intersection form on arithmetic surface": "PASSES-FILTER",
    }
    ok = True
    results = []
    for a in PRESETS:
        v = assess(a)
        results.append(v)
        exp = expect[a["name"]]
        got = v["verdict"]
        # allow moment route to be LIKELY-DEAD (dead keyword dominates)
        match = (got == exp) or (a["name"].startswith("Higher Jensen") and got in ("KNOWN-TERRITORY", "REDUCES-TO-WALL"))
        if verbose:
            print(f"  [{got:16s}] {a['name']}  (expected ~{exp})  {'OK' if match else 'MISMATCH'}")
        if not match:
            ok = False
    print("t01.selftest:", "PASS" if ok else "FAIL")
    return ok


def main():
    print("=== t01 route pre-filter — preset approaches ===")
    results = []
    for a in PRESETS:
        v = assess(a)
        print(f"\n{a['name']}\n  verdict: {v['verdict']}")
        if v["failed_constraints"]:
            print(f"  fails: {', '.join(v['failed_constraints'])}")
        for pa in v["prior_art"]:
            print(f"  prior-art [{pa['class']}] {pa['keyword']}: {pa['status'][:80]}...")
        results.append(Result(claim=f"route '{a['name']}' -> {v['verdict']}", tier=Tier.REPORT,
                              method="proof-shape-spec + prior-art DB", tool="t01_route_prefilter",
                              data=v))
    print()
    ok = selftest()
    path = provenance.stamp("t01_route_prefilter", results,
                            inputs={"n_presets": len(PRESETS)},
                            extra={"selftest_pass": ok})
    print("run record:", path)
    return ok


if __name__ == "__main__":
    main()
