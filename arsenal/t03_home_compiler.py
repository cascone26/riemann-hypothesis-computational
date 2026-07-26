#!/usr/bin/env python3
"""
t03_home_compiler.py  -- Tier S #3.  A SCORER (not a synthesizer) for candidate "homes" for
the RH zeros -- geometric/analytic structures that could host the missing positivity.  Encodes
the Roundtable's 3-clause template (prior / generative / signature-bearing) plus the two hard
gates (reproduce the explicit formula as a Lefschetz trace; carry a definite-signature pairing),
and scores each candidate, ranking them and naming the universal missing clause.

REPORT-tier: it scores structured descriptions against a template drawn from VERIFIED analysis
(rh-roundtable-literature-gate memory + connes-decisive-diagnostic CHECK-8 sec, which read the
full C-C 2006-2026 corpus).  It does NOT build a home (that is the open problem); it makes the
"criteria certify, objects prove" frame operational so we can see which candidate is closest and
exactly what it is missing.

Usage:
    python3 t03_home_compiler.py     # score the known candidate homes + selftest
    from t03_home_compiler import score_home
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core.honesty import Result, Tier
from core import provenance

# The scorecard = Weil's own template (3 clauses) + 2 hard gates.
# Each candidate rates each clause/gate in {2: full, 1: partial, 0: absent}.
CLAUSES = ["prior", "generative", "signature"]
GATES = ["lefschetz_trace", "definite_pairing"]

CLAUSE_TEXT = {
    "prior":      "PRIOR: has an arithmetic-geometry precedent for this kind of positivity.",
    "generative": "GENERATIVE: CONSTRUCTS the object (a home), not just restates a criterion.",
    "signature":  "SIGNATURE-BEARING: carries a definite-signature pairing (Hodge-index), the "
                  "sign source -- NOT a functional-analytic estimate and NOT an identity.",
}
GATE_TEXT = {
    "lefschetz_trace":  "reproduces the explicit formula as a Lefschetz/trace formula.",
    "definite_pairing": "supplies an a-priori definite (Hodge-index) intersection form.",
}

# Candidate homes, rated from the verified corpus read (CHECK-8 + roundtable).
# ratings: {clause/gate: (score 0/1/2, evidence)}
CANDIDATES = {
    "function_field": {
        "label": "Function-field analogue (curve/C over F_q) -- the SOLVED exemplar",
        "prior":            (2, "Weil/Deligne: RH proved via Castelnuovo/Hodge-INDEX, not the trace."),
        "generative":       (2, "The curve and its Jacobian ARE the home; zeros = Frobenius eigenvalues."),
        "signature":        (2, "Hodge-index theorem on the surface CxC gives the definite pairing."),
        "lefschetz_trace":  (2, "Grothendieck-Lefschetz trace formula = the explicit formula."),
        "definite_pairing": (2, "Castelnuovo-Severi inequality = a-priori definite intersection form."),
        "note": "3/3 + both gates. This is WHY it is solved; the template is calibrated to it.",
    },
    "cc_spectral": {
        "label": "Connes-Consani spectral (prolate / scaling Hamiltonian, Weil-positivity)",
        "prior":            (2, "Deep arithmetic-geometry precedent (trace formula, adeles)."),
        "generative":       (2, "Builds a concrete operator/space; reproduces zeros numerically."),
        "signature":        (1, "Positivity = a functional-analytic ESTIMATE (Toeplitz compression), "
                                "all-primes uniformity OPEN. Not a definite pairing. (CHECK-8)"),
        "lefschetz_trace":  (2, "Explicit formula = trace of the scaling action (banked, 1950s+)."),
        "definite_pairing": (0, "No a-priori definite intersection form; the estimate IS the open part."),
        "note": "2/3. Signature clause unmet: estimate, not pairing. The only 2/3 home.",
    },
    "cc_geometric_rr": {
        "label": "Connes-Consani geometric (RR for Spec Z-bar, Jacobian of Spec Z 2026)",
        "prior":            (2, "Arakelov-RR / F1 geometry precedent."),
        "generative":       (2, "Constructs cohomology/theta/Gamma-space, a genuine geometric home."),
        "signature":        (0, "Arakelov-RR IDENTITY tier; 2026 Jacobian builds NO polarization "
                                "(degenerate zero semi-norms). 'positivity remains absent and open'."),
        "lefschetz_trace":  (2, "RR identity plays the trace role."),
        "definite_pairing": (0, "Explicitly no polarization / no definite pairing built (CHECK-8)."),
        "note": "Identity tier: reproduces structure but supplies no sign. Missing clause = signature.",
    },
    "de_branges": {
        "label": "de Branges H(E) space",
        "prior":            (0, "Prior-empty for arithmetic positivity (Conrey-Li: naive version fails)."),
        "generative":       (1, "Constructs a space but A+B has 4 negative eigenvalues."),
        "signature":        (1, "Structure-function inequality is a candidate sign source but unproven."),
        "lefschetz_trace":  (1, "Connects to the functional equation, not a full trace formula."),
        "definite_pairing": (0, "No definite pairing; positivity K(t,t)>0 IS RH, circular."),
        "note": "Prior-empty; not the home.",
    },
    "nyman_beurling": {
        "label": "Nyman-Beurling (closure / approximation criterion)",
        "prior":            (1, "Classical functional-analysis precedent."),
        "generative":       (0, "A CRITERION (density in L^2), not a constructed home -- venue only."),
        "signature":        (0, "No pairing; equivalent restatement of RH."),
        "lefschetz_trace":  (0, "Not trace-shaped."),
        "definite_pairing": (0, "None."),
        "note": "Venue-only: a criterion, not an object. Scores near zero on generative/signature.",
    },
    "borger_f1": {
        "label": "Borger F1 / Lambda-ring geometry",
        "prior":            (2, "Serious F1-geometry precedent."),
        "generative":       (1, "Builds structure but the relevant surface/index is absent."),
        "signature":        (0, "No index / no definite pairing (roundtable: 'Borger no-index')."),
        "lefschetz_trace":  (1, "Partial trace analogue."),
        "definite_pairing": (0, "No index."),
        "note": "No-index: the signature clause is structurally absent.",
    },
}


def score_home(cand: dict):
    clause_scores = {c: cand[c][0] for c in CLAUSES}
    gate_scores = {g: cand[g][0] for g in GATES}
    clauses_full = sum(1 for c in CLAUSES if clause_scores[c] == 2)
    gates_full = sum(1 for g in GATES if gate_scores[g] == 2)
    total = sum(clause_scores.values()) + sum(gate_scores.values())  # out of 10
    missing = [c for c in CLAUSES if clause_scores[c] < 2] + [g for g in GATES if gate_scores[g] < 2]
    return {
        "clauses_full": clauses_full,       # x/3
        "gates_full": gates_full,           # x/2
        "raw_total": total,                 # x/10
        "clause_scores": clause_scores,
        "gate_scores": gate_scores,
        "missing": missing,
        "is_full_home": clauses_full == 3 and gates_full == 2,
    }


def rank():
    scored = []
    for cid, cand in CANDIDATES.items():
        s = score_home(cand)
        scored.append((cid, cand["label"], s))
    scored.sort(key=lambda x: (-x[2]["raw_total"], -x[2]["clauses_full"]))
    return scored


def selftest(verbose=True) -> bool:
    ok = True
    scored = rank()
    # 1. function-field (the solved exemplar) must be the unique full home (3/3 + 2/2)
    fulls = [cid for cid, _, s in scored if s["is_full_home"]]
    if fulls != ["function_field"]:
        ok = False; print("  FAIL: full-home set is not exactly {function_field}:", fulls)
    # 2. cc_spectral must be the top NON-solved candidate and miss exactly the signature/pairing
    nonff = [x for x in scored if x[0] != "function_field"]
    top = nonff[0]
    if top[0] != "cc_spectral":
        ok = False; print("  FAIL: top open candidate is not cc_spectral:", top[0])
    if "signature" not in top[2]["missing"] and "definite_pairing" not in top[2]["missing"]:
        ok = False; print("  FAIL: cc_spectral not flagged as missing signature/pairing")
    # 3. the universal missing clause across ALL open candidates is 'signature' or 'definite_pairing'
    open_cands = [s for cid, _, s in scored if cid != "function_field"]
    if not all(("signature" in s["missing"] or "definite_pairing" in s["missing"]) for s in open_cands):
        ok = False; print("  FAIL: signature/pairing is not the universal gap")
    if verbose:
        print("  rank (raw/10, clauses/3, gates/2):")
        for cid, label, s in scored:
            print(f"    {s['raw_total']:2d}/10  {s['clauses_full']}/3  {s['gates_full']}/2  {label}")
    print("t03.selftest:", "PASS" if ok else "FAIL")
    return ok


def main():
    print("=== t03 Home Compiler — scoring candidate homes for the RH zeros ===")
    print("Template = Weil's 3 clauses (prior/generative/signature) + 2 gates "
          "(Lefschetz trace, definite pairing).\n")
    scored = rank()
    results = []
    for cid, label, s in scored:
        flag = "  <== FULL HOME (solved exemplar)" if s["is_full_home"] else ""
        print(f"[{s['raw_total']:2d}/10]  {label}{flag}")
        print(f"        clauses {s['clause_scores']}  gates {s['gate_scores']}")
        print(f"        missing: {', '.join(s['missing']) or 'none'}")
        print(f"        {CANDIDATES[cid]['note']}")
        results.append(Result(
            claim=f"home '{cid}' scores {s['raw_total']}/10 ({s['clauses_full']}/3 clauses, "
                  f"{s['gates_full']}/2 gates); missing {s['missing']}",
            tier=Tier.REPORT, method="Weil 3-clause scorecard + 2 gates",
            tool="t03_home_compiler", data=s))
    print("\nUniversal gap: every OPEN candidate misses the SIGNATURE clause / definite pairing "
          "-- the a-priori Hodge-index intersection form is the confirmed missing keystone.")
    ok = selftest()
    path = provenance.stamp("t03_home_compiler", results, extra={"selftest_pass": ok})
    print("run record:", path)
    return ok


if __name__ == "__main__":
    main()
