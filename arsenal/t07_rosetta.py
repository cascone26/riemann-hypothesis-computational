#!/usr/bin/env python3
"""
t07_rosetta.py  -- Tier B #7.  The RH-equivalence transport engine.  Holds the family of known
RH-equivalent / RH-implied statements as one linked object and, given a PARTIAL result in one
form, reports its image and relative strength in the others -- so a weak-looking result in one
basis is checked for being near-complete in another.

REPORT-tier knowledge engine.  Edges are curated equivalences/implications from the standard RH
literature and the project's own OBSTRUCTION_MAP; transport annotations say what a partial
result in the source form means downstream.  It does not prove equivalences; it makes the
dictionary queryable so we stop re-deriving that "result X in form A" is "result Y in form B".

Usage:
    python3 t07_rosetta.py
    from t07_rosetta import transport
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import networkx as nx
from core.honesty import Result, Tier
from core import provenance

# equivalent forms of RH (node -> description). All are EQUIVALENT to RH unless noted.
FORMS = {
    "weil_pos":    "Weil positivity: A(g)+B(g) >= P(g) for all admissible bandlimited g",
    "li":          "Li's criterion: lambda_n >= 0 for all n >= 1",
    "debranges":   "de Branges: xi-space structure function is Hermite-Biehler / phase phi'(t)>0",
    "nyman_beur":  "Nyman-Beurling: span of {rho(theta/x)} dense in L^2(0,1)",
    "lagarias":    "Lagarias/Robin: sigma(n) < e^gamma n log log n for n>5040 (=> RH)",
    "hilbert_polya":"Hilbert-Polya: exists self-adjoint H with spectrum = {gamma}",
    "lee_yang":    "Lee-Yang/physics: xi is a limit of Lee-Yang polynomials (reflection+positivity)",
    "jensen":      "Jensen-Polya: all Jensen polynomials J^{d,n}(xi) hyperbolic (uniform in d)",
    "dbn":         "de Bruijn-Newman: Lambda <= 0 (with Rodgers-Tao Lambda>=0 => Lambda=0)",
    "speiser":     "Speiser: zeta'(s) != 0 for 0 < Re s < 1/2",
}

# curated equivalence / implication edges with a transport note.
# (src, dst, relation, transport_note)
LINKS = [
    ("weil_pos", "li", "EQUIV", "Weil form on the Li test functions g_n IS lambda_n >= 0; "
        "a positivity margin for support<=S gives lambda_n>=0 up to n~e^S."),
    ("li", "weil_pos", "EQUIV", "lambda_n = Weil functional at the n-th Li test function."),
    ("weil_pos", "debranges", "EQUIV", "tau=A+B-P PSD <=> reproducing kernel K(t,t)>0 <=> HB."),
    ("debranges", "weil_pos", "EQUIV", "phi'(t)>0 everywhere <=> no Lehmer collision <=> tau PSD."),
    ("weil_pos", "hilbert_polya", "IMPLIES_IF", "a self-adjoint H with the right spectrum makes "
        "tau a Gram matrix => positivity; but no candidate H with correct density exists."),
    ("jensen", "li", "EQUIV_UNIFORM", "uniform-in-d Jensen hyperbolicity <=> all lambda_n>=0; "
        "FIXED-degree Jensen (GORTTW, d<=10^20) is RH-FREE (kernel-local) -- transports to NOTHING "
        "about RH, only to a finite family of Turan inequalities."),
    ("dbn", "jensen", "SAME_WALL", "de Bruijn-Newman heat flow and Jensen hyperbolicity are the "
        "same Lee-Yang wall; Lambda<=0 <=> the flowed xi stays hyperbolic."),
    ("lee_yang", "jensen", "SAME_WALL", "Lee-Yang polynomial limit = Polya-Jensen; reflection is "
        "the functional equation, the MISSING piece is coefficient positivity."),
    ("nyman_beur", "weil_pos", "VENUE_ONLY", "Nyman-Beurling is a criterion (density), not a "
        "constructed home; a partial approximation result transports to a distance bound, not sign."),
    ("lagarias", "li", "SUFFICIENT", "Robin/Lagarias inequality => RH; a partial verification up to "
        "N transports to 'no counterexample below N', a finite check, not the sign mechanism."),
    ("speiser", "weil_pos", "DIFFERENT_AXIS", "Speiser (zeta' zero-free) supports mollifier method: "
        "positive PROPORTION on line (~41.7%); transports to a DENSITY statement, not full positivity; "
        "barrier theta<4/7."),
]


def build():
    G = nx.MultiDiGraph()
    for f, d in FORMS.items():
        G.add_node(f, desc=d)
    for s, d, rel, note in LINKS:
        G.add_edge(s, d, relation=rel, note=note)
    return G


def transport(source_form: str, result_desc: str):
    """Given a partial result stated in `source_form`, list its image/strength in linked forms."""
    G = build()
    if source_form not in FORMS:
        return {"error": f"unknown form {source_form}; known: {list(FORMS)}"}
    out = []
    for _, dst, data in G.out_edges(source_form, data=True):
        out.append({"to": dst, "to_desc": FORMS[dst], "relation": data["relation"],
                    "transport": data["note"]})
    return {"source": source_form, "source_desc": FORMS[source_form],
            "result": result_desc, "images": out}


def selftest(verbose=True) -> bool:
    ok = True
    # 1. Jensen fixed-degree must transport to 'RH-FREE / nothing about RH' (the literature-gate lesson)
    t = transport("jensen", "hyperbolic for all d<=10^20 (GORTTW)")
    jl = [im for im in t["images"] if im["to"] == "li"]
    if not jl or "RH-FREE" not in jl[0]["transport"]:
        ok = False; print("  FAIL: fixed-degree Jensen did not transport to RH-FREE")
    # 2. Weil positivity <-> Li must be marked EQUIV both ways
    a = transport("weil_pos", "margin>0 for support<=S")
    if not any(im["to"] == "li" and im["relation"] == "EQUIV" for im in a["images"]):
        ok = False; print("  FAIL: weil<->li not EQUIV")
    # 3. Speiser transports to a DIFFERENT_AXIS density statement, not full positivity
    s = transport("speiser", "positive proportion 41.7% on line")
    if not any(im["relation"] == "DIFFERENT_AXIS" for im in s["images"]):
        ok = False; print("  FAIL: speiser not flagged different-axis")
    if verbose:
        print("  transport('jensen', fixed-degree result):")
        for im in t["images"]:
            print(f"    -> {im['to']} [{im['relation']}]: {im['transport'][:90]}...")
    print("t07.selftest:", "PASS" if ok else "FAIL")
    return ok


def main():
    print("=== t07 Rosetta: RH-equivalence transport engine ===\n")
    G = build()
    print(f"forms: {G.number_of_nodes()}, links: {G.number_of_edges()}\n")
    demos = [
        ("jensen", "Jensen polynomials hyperbolic for all fixed d <= 10^20 (GORTTW)"),
        ("speiser", "positive proportion 41.7% of zeros on the line (Levinson-Conrey)"),
        ("weil_pos", "positivity margin > 0 for test-function support <= C (first theorem target)"),
    ]
    results = []
    for form, desc in demos:
        t = transport(form, desc)
        print(f"RESULT in [{form}]: {desc}")
        for im in t["images"]:
            print(f"   => [{im['to']}] ({im['relation']})")
            print(f"      {im['transport']}")
        print()
        results.append(Result(
            claim=f"transported '{desc[:50]}...' from {form} to {len(t['images'])} linked forms",
            tier=Tier.REPORT, method="curated RH-equivalence graph", tool="t07_rosetta", data=t))
    ok = selftest()
    path = provenance.stamp("t07_rosetta", results, extra={"selftest_pass": ok})
    print("run record:", path)
    return ok


if __name__ == "__main__":
    main()
