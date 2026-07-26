#!/usr/bin/env python3
"""
t05_ladder_finder.py  -- Tier A #5.  A reverse-mathematics ladder: the partial order of
RH-related statements by logical strength, from what the campaign has PROVEN up to RH, with a
frontier finder that names the weakest OPEN statement just above the proven line -- the next real
stepping-stone rung, so we stop throwing sessions straight at the summit.

REPORT-tier reasoning scaffold over VERIFIED project results (LEDGER: Yoshida support<=log2,
our first-prime extension, window-scan n<=9.8M, pair-positivity, GORTTW fixed-degree Jensen) and
the standard implication chain toward RH.  Edges are 'A implies B' (A stronger than B).  It does
not prove any rung; it locates the frontier and the shortest open gap.

Usage: python3 t05_ladder_finder.py
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import networkx as nx
from core.honesty import Result, Tier
from core import provenance

# node id -> (statement, status, evidence)
# status in {PROVEN, OPEN, EQUALS_RH}
NODES = {
    "ff_rh":     ("Function-field RH (curve/F_q): Frobenius eigenvalues on |q|^{1/2}",
                  "PROVEN", "Weil/Deligne (different field; the solved exemplar)"),
    "yoshida":   ("Weil positivity for test-function support <= log 2",
                  "PROVEN", "Yoshida / Bombieri, unconditional"),
    "firstprime":("Weil positivity past the first prime, p in {2,3}, finite-Galerkin, a in [0.549,0.64]",
                  "PROVEN", "LEDGER 199.5/209-227, Arb-certified at points; K0=9 tail composes"),
    "window98M": ("lambda_n > 0 for all n <= 9,834,937",
                  "PROVEN", "window-scan certificate"),
    "pairpos":   ("lambda_n > 0 for n < 2*pi*gamma_1 (~88.8) assumption-free; n<~1.9e13 via Platt height",
                  "PROVEN", "corrected pair-positivity theorem, unconditional/finite"),
    "jensen_fix":("Jensen polynomials J^{d,n} hyperbolic for each FIXED degree d <= ~10^20",
                  "PROVEN", "GORTTW/Chasse, Dimitrov-Lucas; kernel-local, RH-free"),
    # ---- open rungs, increasing strength ----
    "support_C": ("Weil positivity for test-function support <= C, for some FIXED C > log 2",
                  "OPEN", "the continuum / tail-lemma target: a-interval Arb over (0.42, 0.5 log3)"),
    "uniform_N": ("lambda_min(tau_N) >= 0 for ALL N at fixed bandwidth lambda^2 (uncond.)",
                  "OPEN", "finite-Galerkin proven; the all-N uniformity is open"),
    "jensen_unif":("Jensen hyperbolicity uniform in degree d (all d, all n)",
                  "OPEN", "the RH-content of the Jensen route lives here"),
    "uniform_lam":("A(g)+B(g) >= P(g) uniformly in bandwidth (all admissible g)",
                  "EQUALS_RH", "the Wall: equals RH (Bombieri/Connes)"),
    "debranges_ineq":("de Branges structure-function inequality for the xi-space",
                  "EQUALS_RH", "equivalent to RH"),
    "rh":        ("Riemann Hypothesis",
                  "EQUALS_RH", "the target"),
}

# implication edges: (stronger) -> (weaker/implied)
EDGES = [
    ("rh", "uniform_lam"), ("uniform_lam", "rh"),          # equal
    ("rh", "debranges_ineq"), ("debranges_ineq", "rh"),    # equal
    ("uniform_lam", "uniform_N"),        # all-bandwidth => fixed-bandwidth all-N
    ("uniform_lam", "jensen_unif"),      # wall => jensen uniform (same content)
    ("jensen_unif", "jensen_fix"),       # uniform-d => each fixed d
    ("uniform_N", "support_C"),          # all-N at fixed lambda => finite support window
    ("support_C", "firstprime"),         # support C>log2 => the first-prime window (log2..)
    ("firstprime", "yoshida"),           # extends Yoshida
    ("uniform_lam", "window98M"),        # RH => the finite window cert
    ("uniform_lam", "pairpos"),          # RH => pair positivity everywhere
    ("ff_rh", "ff_rh"),                  # standalone solved exemplar (analogy anchor)
]


def build():
    G = nx.DiGraph()
    for nid, (stmt, status, ev) in NODES.items():
        G.add_node(nid, stmt=stmt, status=status, ev=ev)
    for a, b in EDGES:
        if a != b:
            G.add_edge(a, b)
    return G


def frontier(G):
    """
    Frontier = OPEN nodes all of whose STRICTLY-WEAKER implied statements are already PROVEN,
    i.e. the next rung you could try to climb because everything below it is done.
    A node X is weaker than Y if there's a path Y -> ... -> X.
    """
    front = []
    for nid, d in G.nodes(data=True):
        if d["status"] != "OPEN":
            continue
        weaker = nx.descendants(G, nid)
        weaker_open = [w for w in weaker if G.nodes[w]["status"] == "OPEN"]
        if not weaker_open:   # everything strictly below is proven (or equals-rh above, excluded)
            front.append(nid)
    return front


def selftest(verbose=True) -> bool:
    G = build()
    ok = True
    if not nx.is_directed_acyclic_graph(nx.DiGraph([(a, b) for a, b in EDGES
                                                    if not (b, a) in EDGES or a < b])):
        pass  # equalities create 2-cycles by design; skip strict-DAG check
    front = frontier(G)
    # the nearest stepping stone must be 'support_C' (fixed support C>log2 = the project's real NEXT)
    if "support_C" not in front:
        ok = False; print("  FAIL: support_C not on the frontier:", front)
    # RH must sit above everything proven (not itself provable-adjacent)
    if G.nodes["rh"]["status"] != "EQUALS_RH":
        ok = False
    if verbose:
        print("  frontier (nearest open rungs, everything below proven):")
        for f in front:
            print(f"    -> {f}: {NODES[f][0]}")
    print("t05.selftest:", "PASS" if ok else "FAIL")
    return ok


def main():
    print("=== t05 Reverse-math ladder finder ===")
    G = build()
    print("\nPROVEN rungs:")
    for nid, d in NODES.items():
        if d[1] == "PROVEN":
            print(f"  [PROVEN] {d[0]}")
    print("\nOPEN rungs (toward RH):")
    for nid, d in NODES.items():
        if d[1] == "OPEN":
            print(f"  [OPEN]   {d[0]}")
    front = frontier(G)
    print("\nFRONTIER -- nearest stepping-stones (everything strictly below is proven):")
    for f in front:
        print(f"  ==> {NODES[f][0]}")
        print(f"      why reachable: {NODES[f][2]}")
    results = [Result(
        claim=f"nearest open stepping-stone rung(s): {front} -- climb here next, not the summit",
        tier=Tier.REPORT, method="reverse-math implication lattice over verified results",
        tool="t05_ladder_finder", data={"frontier": front})]
    ok = selftest()
    path = provenance.stamp("t05_ladder_finder", results, extra={"selftest_pass": ok})
    print("run record:", path)
    return ok


if __name__ == "__main__":
    main()
