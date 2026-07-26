#!/usr/bin/env python3
"""
t02_wall_atlas.py  -- Tier S #2.  A machine-readable graph of every mapped RH route, its
precise obstruction, and the reductions between them, with an equivalence detector that finds
routes dying for the SAME reason (the Wall Invariant, made mechanical).

REPORT-tier structural tool over already-verified data (OBSTRUCTION_MAP.md + RH_WALL_ANATOMY.md
sec 3).  It does not discover new obstructions; it makes the campaign's own conclusion --
"every classical route reduces to the missing sign-carrying bridge" -- checkable and queryable,
so no future session re-treads and the exact common root is explicit.

Usage:
    python3 t02_wall_atlas.py            # build graph, print equivalence classes + common root
    from t02_wall_atlas import build_atlas, equivalence_classes
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import networkx as nx
from core.honesty import Result, Tier
from core import provenance

# route id -> (label, obstruction node, status)
ROUTES = {
    "moment":        ("Moment / Hankel (Track 1)",        "obs_logconcave",   "DEAD"),
    "beurling":      ("Beurling-Selberg majorant",        "obs_pointwise",    "WALL"),
    "large_sieve":   ("Large sieve global bound",         "obs_offline_only", "WALL"),
    "debranges":     ("de Branges H(E) space",            "obs_neg_eigs",     "WALL"),
    "operator_tau":  ("Operator / tau spectral (Track 2)","obs_is_weil",      "WALL"),
    "pair_pos":      ("Pair positivity + verif height",   "obs_finite_height","PARTIAL"),
    "speiser":       ("Speiser / Levinson mollifier",     "obs_theta_47",     "PARTIAL"),
    "jensen":        ("Jensen / Polya-Laguerre",          "obs_fixed_degree", "KNOWN"),
    "debruijn":      ("de Bruijn-Newman detector",        "obs_detector_only","WALL"),
    "crystalline":   ("Crystalline / Fourier-quasicrystal","obs_leeyang",     "WALL"),
}

# obstruction node -> (human text, reduces_to)  ; reduces_to = None means terminal/distinct
OBSTRUCTIONS = {
    "obs_logconcave":   ("lambda_n log-concave, not a moment sequence (det<0 at N=2)", "WALL_INVARIANT"),
    "obs_pointwise":    ("prime comb exceeds arch kernel at every prime>=5; no pointwise majorant", "WALL_INVARIANT"),
    "obs_offline_only": ("only bounds off-line zeros (zero-density); never forces them on", "WALL_INVARIANT"),
    "obs_neg_eigs":     ("A+B has 4 negative eigenvalues; only tau=A+B-P is PSD", "WALL_INVARIANT"),
    "obs_is_weil":      ("every tau certificate reduces to Weil; tau>=0 IS RH", "WALL_INVARIANT"),
    "obs_fixed_degree": ("fixed-degree Jensen ladder is kernel-local & RH-free; RH only in uniformity in d", "WALL_INVARIANT"),
    "obs_detector_only":("g_lambda is zero-side only; no Euler product => no prime handle", "WALL_INVARIANT"),
    "obs_leeyang":      ("FQ classification passes through Lee-Yang polys = Jensen wall", "WALL_INVARIANT"),
    # the two that do NOT reduce to the wall (genuinely different axes) -- honest distinction
    "obs_finite_height":("unconditional but finite: lambda_n>0 for n<~2pi*T; not the wall, just finite", None),
    "obs_theta_47":     ("shifted-convolution barrier theta<4/7; different axis, not the positivity wall", None),
}

WALL_INVARIANT = ("WALL_INVARIANT",
    "Every classical route reduces to a required INEQUALITY between a prime-local object and a "
    "zero-global object. The only prime<->zero identity is the explicit formula, and it is a "
    "sign-blind EQUALITY. The classical toolkit contains NO second, sign-carrying bridge.")


def build_atlas() -> nx.DiGraph:
    G = nx.DiGraph()
    G.add_node(WALL_INVARIANT[0], kind="invariant", text=WALL_INVARIANT[1])
    for rid, (label, obs, status) in ROUTES.items():
        G.add_node(rid, kind="route", label=label, status=status)
        G.add_edge(rid, obs, rel="obstructed_by")
    for oid, (text, reduces) in OBSTRUCTIONS.items():
        G.add_node(oid, kind="obstruction", text=text)
        if reduces:
            G.add_edge(oid, reduces, rel="reduces_to")
    return G


def equivalence_classes(G: nx.DiGraph):
    """Routes are equivalent iff they reduce to the same terminal node (obstruction or invariant)."""
    terminal_of = {}
    for rid in ROUTES:
        # follow reduces_to edges to a terminal
        obs = ROUTES[rid][1]
        node = obs
        seen = set()
        while True:
            succ = [v for _, v, d in G.out_edges(node, data=True) if d.get("rel") == "reduces_to"]
            if not succ or node in seen:
                break
            seen.add(node)
            node = succ[0]
        terminal_of[rid] = node
    classes = {}
    for rid, term in terminal_of.items():
        classes.setdefault(term, []).append(rid)
    return classes, terminal_of


def selftest(verbose=True) -> bool:
    G = build_atlas()
    ok = True
    # 1. graph connected as a DAG rooted logically at WALL_INVARIANT for the wall routes
    if not nx.is_directed_acyclic_graph(G):
        ok = False; print("  FAIL: atlas has a cycle")
    classes, term = equivalence_classes(G)
    wall_routes = classes.get("WALL_INVARIANT", [])
    # 2. the five classical wall routes (2-5 + jensen/crystalline/debruijn) collapse to the invariant
    expected_wall = {"moment", "beurling", "large_sieve", "debranges", "operator_tau",
                     "jensen", "debruijn", "crystalline"}
    if not expected_wall.issubset(set(wall_routes)):
        ok = False; print("  FAIL: not all wall routes reduce to WALL_INVARIANT:",
                          expected_wall - set(wall_routes))
    # 3. pair_pos and speiser are honestly NOT in the wall class (distinct axes)
    if "pair_pos" in wall_routes or "speiser" in wall_routes:
        ok = False; print("  FAIL: a distinct-axis route wrongly folded into the wall")
    if verbose:
        print(f"  routes: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
        print(f"  equivalence classes ({len(classes)}):")
        for term, members in classes.items():
            print(f"    -> {term}: {', '.join(members)}")
    print("t02.selftest:", "PASS" if ok else "FAIL")
    return ok


def main():
    G = build_atlas()
    classes, term = equivalence_classes(G)
    print("=== t02 Wall Atlas ===")
    print(f"nodes={G.number_of_nodes()} edges={G.number_of_edges()}\n")
    print("WALL INVARIANT (common root):")
    print("  " + WALL_INVARIANT[1] + "\n")
    print("Equivalence classes (routes dying for the same reason):")
    for t, members in sorted(classes.items(), key=lambda kv: -len(kv[1])):
        labels = [ROUTES[m][0] for m in members]
        print(f"  [{t}]  ({len(members)} routes)")
        for m in members:
            print(f"      - {ROUTES[m][0]}  ({ROUTES[m][2]})")
    print()
    ok = selftest()
    results = [Result(
        claim=f"{len(classes.get('WALL_INVARIANT',[]))} classical routes reduce to the single "
              f"Wall Invariant (no sign-carrying bridge); 2 routes (pair-pos, mollifier) are "
              f"honestly distinct finite axes",
        tier=Tier.REPORT, method="reduction graph over OBSTRUCTION_MAP + RH_WALL_ANATOMY",
        tool="t02_wall_atlas",
        data={"classes": {k: v for k, v in classes.items()}})]
    path = provenance.stamp("t02_wall_atlas", results, extra={"selftest_pass": ok})
    print("run record:", path)
    return ok


if __name__ == "__main__":
    main()
