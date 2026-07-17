#!/usr/bin/env python3
"""
community_detect.py — PRISM EXTREME item #6: Complex-networks community detection.

Builds the eigenvector-correlation graph on the near-zero eigenspace of tau_N
across multiple (N, lsq) pairs, then runs community detection.

Graph construction:
  Nodes  = (N, k) pairs — the k-th eigenvector of tau_N at a given lsq.
  Edges  = |corr(v_k^N, v_l^M)|  between eigenvector component vectors.
  Weight threshold  θ = 0.3 (only strong correlations).

Community detection: greedy modularity maximization (scipy/networkx-free version
using the simple label-propagation scheme).

Physical interpretation:
  - A single community = all near-zero eigenvectors are the same "mode" (stable ground state)
  - Multiple communities = distinct mode families (RSB or multiple branches of the
    zero-manifold)
  - Community structure that clusters by N (not by index k) → N-evolution is
    continuous (adiabatic) → no phase transition between N values
  - Community structure that clusters by k (not by N) → mode families are stable
    across N (persistent spectral structures)
"""

import math
import numpy as np
import scipy.linalg
from scipy.special import digamma as _digamma
import json
import time

TWO_PI = 2*math.pi
LOGPI = math.log(math.pi)

def omega_vec(r):
    return np.real(_digamma(0.25+0.5j*r)) - LOGPI

def build_sieve(lsq):
    sieve = np.ones(lsq+1, dtype=bool)
    sieve[0] = sieve[1] = False
    for i in range(2, int(lsq**0.5)+1):
        if sieve[i]: sieve[i*i::i] = False
    lkpp, wts = [], []
    for p in range(2, lsq+1):
        if not sieve[p]: continue
        pk = p
        while pk <= lsq:
            lkpp.append(math.log(pk))
            wts.append(math.log(p)/math.sqrt(pk))
            pk *= p
    return np.array(lkpp), np.array(wts)

def build_tau(N, lsq, G=6000, K_BATCH=2000):
    L = 2*math.log(math.sqrt(lsq))
    ns = np.arange(-N, N+1)
    D = len(ns)
    R = TWO_PI*N/L + 60.0
    grid = np.linspace(-R, R, G)
    h = grid[1]-grid[0]
    w = np.full(G, h); w[0]=w[-1]=0.5*h
    Om = omega_vec(grid)
    c = TWO_PI*ns[:,None]/L
    d = c - grid[None,:]
    sm = np.abs(d)<1e-12
    E = np.where(sm, np.sqrt(L)+0j,
                 (1/np.sqrt(L))*np.exp(0.5j*d*L)*2*np.sin(0.5*d*L)/d)
    arch = (E*(w*Om) @ E.conj().T).real / TWO_PI
    sh = np.sinh(L/4)**2
    ni = ns[:,None]; nj = ns[None,:]
    num = L**2 - 16*np.pi**2*ni*nj
    den = (L**2+16*np.pi**2*ni**2)*(L**2+16*np.pi**2*nj**2)
    B = 32*L*sh*num/den
    lkpp, wts_p = build_sieve(lsq)
    K = len(lkpp)
    S = np.zeros(D); S_diag = np.zeros(D)
    for k0 in range(0, K, K_BATCH):
        k1 = min(k0+K_BATCH, K)
        ph = TWO_PI*np.outer(ns, lkpp[k0:k1])/L
        S += np.sin(ph) @ wts_p[k0:k1]
        S_diag += np.cos(ph) @ (2*wts_p[k0:k1]*(1-lkpp[k0:k1]/L))
    Si = S[:,None]; Sj = S[None,:]
    dni = ns[:,None] - ns[None,:]
    with np.errstate(divide='ignore', invalid='ignore'):
        P = np.where(dni != 0, (Sj-Si)/(math.pi*dni), 0.0)
    np.fill_diagonal(P, S_diag)
    tau = arch+B-P
    return (tau+tau.T)/2

def label_propagation(adj, n_iters=100, seed=42):
    """Simple label propagation community detection."""
    rng = np.random.default_rng(seed)
    n = adj.shape[0]
    labels = np.arange(n)
    for _ in range(n_iters):
        order = rng.permutation(n)
        changed = False
        for i in order:
            neighbors = np.where(adj[i] > 0)[0]
            if len(neighbors) == 0:
                continue
            # Most frequent label among neighbors (weighted)
            neighbor_labels = labels[neighbors]
            neighbor_weights = adj[i, neighbors]
            label_scores = {}
            for l, w in zip(neighbor_labels, neighbor_weights):
                label_scores[l] = label_scores.get(l, 0) + w
            best_label = max(label_scores, key=label_scores.get)
            if best_label != labels[i]:
                labels[i] = best_label
                changed = True
        if not changed:
            break
    # Relabel to 0..n_communities-1
    unique = sorted(set(labels))
    remap = {l: i for i, l in enumerate(unique)}
    return np.array([remap[l] for l in labels])

def modularity(adj, communities):
    """Compute Newman-Girvan modularity Q."""
    n = adj.shape[0]
    m = adj.sum() / 2
    if m == 0:
        return 0.0
    k = adj.sum(axis=1)  # degree vector
    Q = 0.0
    for i in range(n):
        for j in range(i, n):
            if communities[i] == communities[j]:
                e_ij = adj[i,j] if i != j else 0
                expected = k[i]*k[j] / (2*m)
                Q += (e_ij - expected) * (2 if i != j else 1)
    return Q / (2*m)

if __name__ == "__main__":
    print("Community Detection on Eigenvector-Correlation Network")
    print("=========================================================")
    print()

    # Compute near-zero eigenvectors at multiple (N, lsq) pairs
    # Focus on two lsq values: 50 (near-resonant) and 14 (non-resonant)
    configs = []
    for lsq in [14, 50]:
        for N in [6, 8, 10, 12, 16, 20]:
            configs.append((N, lsq))

    print(f"Building {len(configs)} operator instances...")
    print()

    nodes = []   # (N, lsq, k, vec) — each near-zero eigenvector is a node
    node_meta = []

    MAX_D = 2*20+1  # 41 — pad all eigenvectors to this size

    for N, lsq in configs:
        t0 = time.time()
        tau = build_tau(N, lsq)
        D = 2*N+1
        evals, evecs = scipy.linalg.eigh(tau)
        shift = max(0, -evals[0] + 1e-10)
        evals_s = evals + shift

        # Find spectral gap to identify "near-zero" block
        if D > 2:
            gaps = np.diff(evals_s)
            gap_idx = np.argmax(gaps[:D//2])
            nz_idx = np.arange(gap_idx+1)
        else:
            nz_idx = np.array([0])

        for k in nz_idx:
            v = evecs[:, k]
            # Pad/center to MAX_D
            v_pad = np.zeros(MAX_D)
            offset = (MAX_D - D) // 2
            v_pad[offset:offset+D] = v
            nodes.append(v_pad)
            node_meta.append({
                "N": N, "lsq": lsq, "k": int(k),
                "ev": float(evals[k]), "shift": float(shift)
            })
        dt = time.time()-t0
        print(f"  N={N}, lsq={lsq}: D={D}, n_nz={len(nz_idx)}, "
              f"min_ev={evals[0]:.4e}  ({dt:.1f}s)")

    print()
    n_nodes = len(nodes)
    print(f"Total nodes: {n_nodes}")
    if n_nodes == 0:
        print("No near-zero eigenvectors found.")
        import sys; sys.exit(1)

    # Build adjacency matrix: |correlation| between all pairs
    X = np.array(nodes)  # n_nodes × MAX_D
    # Normalize each row
    norms = np.linalg.norm(X, axis=1, keepdims=True)
    norms = np.where(norms < 1e-12, 1.0, norms)
    X_norm = X / norms

    # Correlation matrix (Gram matrix in normalized basis)
    C = np.abs(X_norm @ X_norm.T)  # n_nodes × n_nodes

    # Threshold: only keep edges where |corr| > 0.3
    THETA = 0.3
    adj = np.where(C > THETA, C, 0.0)
    np.fill_diagonal(adj, 0.0)  # no self-loops

    n_edges = (adj > 0).sum() // 2
    print(f"Graph: {n_nodes} nodes, {n_edges} edges (|corr|>{THETA})")

    # Connectivity analysis
    # Connected components via BFS
    visited = np.zeros(n_nodes, dtype=bool)
    components = []
    for start in range(n_nodes):
        if visited[start]:
            continue
        comp = []
        queue = [start]
        while queue:
            node = queue.pop(0)
            if visited[node]:
                continue
            visited[node] = True
            comp.append(node)
            neighbors = np.where(adj[node] > 0)[0]
            queue.extend(neighbors)
        components.append(comp)

    print(f"Connected components: {len(components)}")
    for i, comp in enumerate(sorted(components, key=len, reverse=True)[:10]):
        ns_in = sorted(set(node_meta[j]['N'] for j in comp))
        lsqs_in = sorted(set(node_meta[j]['lsq'] for j in comp))
        print(f"  Component {i+1}: {len(comp)} nodes, N∈{ns_in}, lsq∈{lsqs_in}")

    print()
    print("Running label propagation community detection...")
    communities = label_propagation(adj)
    n_communities = len(set(communities))
    print(f"Communities found: {n_communities}")

    # Modularity
    Q = modularity(adj, communities)
    print(f"Modularity Q = {Q:.4f}  "
          f"({'strong communities Q>0.3' if Q > 0.3 else 'weak communities Q<0.3'})")

    print()
    print("Community breakdown:")
    for c_id in range(n_communities):
        members = [i for i in range(n_nodes) if communities[i] == c_id]
        if not members:
            continue
        N_vals = [node_meta[i]['N'] for i in members]
        lsq_vals = [node_meta[i]['lsq'] for i in members]
        k_vals = [node_meta[i]['k'] for i in members]
        print(f"  Community {c_id}: {len(members)} members")
        print(f"    N ∈ {sorted(set(N_vals))}")
        print(f"    lsq ∈ {sorted(set(lsq_vals))}")
        print(f"    k ∈ {sorted(set(k_vals))}")

        # Homogeneity analysis: does this community cluster by N or by k?
        N_entropy = 0.0
        from collections import Counter
        for val, count in Counter(N_vals).items():
            p = count / len(members)
            N_entropy -= p * np.log(p)
        k_entropy = 0.0
        for val, count in Counter(k_vals).items():
            p = count / len(members)
            k_entropy -= p * np.log(p)
        print(f"    Entropy: N_entropy={N_entropy:.3f}, k_entropy={k_entropy:.3f}  "
              f"({'clusters by k (mode)' if k_entropy < N_entropy else 'clusters by N (depth)'})")

    print()
    # Check if lsq=14 and lsq=50 form separate communities
    lsq14_nodes = [i for i in range(n_nodes) if node_meta[i]['lsq'] == 14]
    lsq50_nodes = [i for i in range(n_nodes) if node_meta[i]['lsq'] == 50]
    lsq14_comms = set(communities[i] for i in lsq14_nodes)
    lsq50_comms = set(communities[i] for i in lsq50_nodes)
    overlap = lsq14_comms & lsq50_comms
    print(f"lsq=14 communities: {sorted(lsq14_comms)}")
    print(f"lsq=50 communities: {sorted(lsq50_comms)}")
    print(f"Shared communities: {sorted(overlap)}")
    if overlap:
        print("  → lsq=14 and lsq=50 eigenvectors share communities: CORRELATED across lsq")
    else:
        print("  → lsq=14 and lsq=50 eigenvectors form SEPARATE communities: DECOUPLED")

    print()
    print("=== PHYSICAL INTERPRETATION ===")
    print(f"  n_communities = {n_communities}")
    if n_communities == 1:
        print("  Single community: all near-zero modes are the same family.")
        print("  → Simple ground state structure (consistent with unique RH zero line).")
    elif n_communities <= 3:
        print(f"  {n_communities} distinct communities: a small number of mode families.")
        print("  → Structured near-zero manifold (possibly k=0,±1,...  angular momentum sectors).")
    else:
        print(f"  {n_communities} communities: rich multi-branch structure.")
        print("  → RSB-like near-zero manifold (many distinct mode families → complex landscape).")

    if Q > 0.3:
        print(f"  Strong modularity (Q={Q:.3f}): community structure is statistically robust.")
    else:
        print(f"  Weak modularity (Q={Q:.3f}): community structure may be noise.")

    result = {
        "n_nodes": n_nodes,
        "n_edges": int(n_edges),
        "n_components": len(components),
        "n_communities": n_communities,
        "modularity_Q": Q,
        "community_assignments": communities.tolist(),
        "node_meta": node_meta,
        "theta": THETA
    }
    out = "/Users/scones/projects/riemann/results/community_detect.json"
    with open(out, "w") as fh:
        json.dump(result, fh, indent=2)
    print(f"\nSaved → {out}")
