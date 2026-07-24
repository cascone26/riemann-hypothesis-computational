"""
Group C — real-world network spectra (SNAP email-Eu-core, ~1000 nodes,
25571 directed edges, real institutional email communication network).
Tests: (1) real adjacency/Laplacian (undirected, real symmetric -> expect
GOE if RMT-universal at all, given strong community structure might break
universality), (2) twisted/magnetic complex-Hermitian version (expect GUE,
extending the Group A synthetic-graph finding to a genuine empirical
network for the first time in this project).
"""
import os
import sys
import json
import numpy as np
import networkx as nx

sys.path.insert(0, os.path.dirname(__file__))
import lib
from group_a_twisted import twisted_adjacency_eigs

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
                            "results", "pattern_hunt")


def load_graph(path="/tmp/email-eu-core.txt"):
    G = nx.Graph()
    with open(path) as f:
        for line in f:
            a, b = line.split()
            if a != b:
                G.add_edge(int(a), int(b))
    return G


def main():
    G = load_graph()
    print(f"email-Eu-core: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    candidates = {}
    A = nx.adjacency_matrix(G).toarray().astype(float)
    L = nx.laplacian_matrix(G).toarray().astype(float)
    candidates["realnet_email_adjacency"] = np.linalg.eigvalsh(A)
    candidates["realnet_email_laplacian"] = np.linalg.eigvalsh(L)
    candidates["realnet_email_twisted"] = twisted_adjacency_eigs(G, seed=100)

    results = {}
    for name, seq in candidates.items():
        seq = np.asarray(seq, dtype=float)
        seq = seq[np.isfinite(seq)]
        r = lib.mean_ratio(seq)
        results[name] = {"n": len(seq), "mean_ratio": r}
        print(lib.summarize({"name": name, "n": len(seq), "mean_ratio": r}))

    outpath = os.path.join(RESULTS_DIR, "group_c_realnet_ratios.json")
    with open(outpath, "w") as f:
        json.dump(results, f, indent=2, default=str)

    raw_out = {name: candidates[name].tolist() for name in candidates}
    with open(os.path.join(RESULTS_DIR, "group_c_realnet_raw.json"), "w") as f:
        json.dump(raw_out, f)
    print(f"\nSaved to {outpath}")
    return results


if __name__ == "__main__":
    main()
