"""
Twisted (magnetic) graph spectra — the standard RMT recipe to get GUE
(not GOE) out of a graph: break time-reversal symmetry by putting a random
U(1) phase on each directed edge of a complex Hermitian adjacency matrix,
H_ij = exp(i*theta_ij), H_ji = conj(H_ij) (theta_ij ~ Uniform(0, 2pi) i.i.d.
per undirected edge). This is the discrete analogue of a magnetic
Laplacian / Harper-like operator and is documented to show GUE statistics
once flux is generic (unlike the real adjacency matrix, which is GOE).
Not found already tested against zeta zeros in our own project history —
genuinely a new comparison for us to run.
"""
import os
import sys
import json
import numpy as np
import networkx as nx

sys.path.insert(0, os.path.dirname(__file__))
import lib

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
                            "results", "pattern_hunt")


def twisted_adjacency_eigs(G, seed=0):
    rng = np.random.default_rng(seed)
    n = G.number_of_nodes()
    idx = {node: i for i, node in enumerate(G.nodes())}
    H = np.zeros((n, n), dtype=complex)
    for u, v in G.edges():
        theta = rng.uniform(0, 2 * np.pi)
        i, j = idx[u], idx[v]
        H[i, j] = np.exp(1j * theta)
        H[j, i] = np.exp(-1j * theta)
    eigs = np.linalg.eigvalsh(H)
    return eigs


def main():
    candidates = {}
    print("Building twisted (magnetic) graph spectra...")
    for d in [3, 4, 6, 10]:
        G = nx.random_regular_graph(d, 1200, seed=1)
        candidates[f"twisted_regular_d{d}_N1200"] = twisted_adjacency_eigs(G, seed=10 + d)

    for p in [0.01, 0.05, 0.2]:
        G = nx.erdos_renyi_graph(1200, p, seed=2)
        candidates[f"twisted_erdos_renyi_p{p}_N1200"] = twisted_adjacency_eigs(G, seed=20)

    G = nx.watts_strogatz_graph(1200, 6, 0.3, seed=3)
    candidates["twisted_watts_strogatz_N1200"] = twisted_adjacency_eigs(G, seed=30)

    G = nx.barabasi_albert_graph(1200, 4, seed=4)
    candidates["twisted_barabasi_albert_N1200"] = twisted_adjacency_eigs(G, seed=40)

    results = {}
    for name, seq in candidates.items():
        seq = np.asarray(seq, dtype=float)
        seq = seq[np.isfinite(seq)]
        r = lib.mean_ratio(seq)
        results[name] = {"n_raw": len(seq), "mean_ratio": r}
        print(lib.summarize({"name": name, "n": len(seq), "mean_ratio": r}))

    outpath = os.path.join(RESULTS_DIR, "group_a_twisted_ratios.json")
    with open(outpath, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved to {outpath}")

    raw_out = {name: candidates[name].tolist() for name in candidates}
    with open(os.path.join(RESULTS_DIR, "group_a_twisted_raw.json"), "w") as f:
        json.dump(raw_out, f)

    return results


if __name__ == "__main__":
    main()
