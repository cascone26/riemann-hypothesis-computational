"""
Group B/C candidates computable locally without network access:
  - Circular quantum billiard spectrum (Bessel zeros squared) — INTEGRABLE,
    Berry-Tabor predicts Poisson. Clean textbook negative-control-by-physics.
  - Airy zeros (Bessel-family special function zeros) — untested wildcard.
  - Fibonacci quasicrystal atom positions (cut-and-project) — deterministic
    two-gap structure, expected NOT GUE (hyperuniform, not level-repulsive
    in the RMT sense).
"""
import os
import sys
import math
import json
import numpy as np
from scipy import special

sys.path.insert(0, os.path.dirname(__file__))
import lib

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
                            "results", "pattern_hunt")


def circular_billiard_spectrum(m_max=200, k_max=200):
    """
    Dirichlet Laplacian eigenvalues on the unit disk: E_{m,k} = j_{m,k}^2
    where j_{m,k} is the k-th positive zero of Bessel function J_m.
    The circular billiard is INTEGRABLE -> Berry-Tabor conjecture predicts
    Poisson statistics (no level repulsion), unlike chaotic billiards
    (Bunimovich stadium, Sinai) which show GOE/GUE.
    """
    eigs = []
    for m in range(m_max):
        zeros = special.jn_zeros(m, k_max)
        eigs.extend((zeros ** 2).tolist())
        # degeneracy: m>0 modes are doubly degenerate (cos/sin), but for a
        # spectral-statistics test we use DISTINCT eigenvalues (a physical
        # measurement would resolve or not resolve degeneracy depending on
        # symmetry-breaking; report both variants for honesty)
    return np.array(sorted(eigs))


def airy_zeros(n=5000):
    """Zeros of the Airy function Ai(x) — untested wildcard special-function
    zero sequence (unlike Bessel, no direct billiard/physical interpretation
    used here; purely a numerical curiosity check)."""
    m = np.arange(1, n + 1)
    zeros = special.ai_zeros(n)[0]  # returns (a_ai, ap_ai, ai_val, aip_val)
    return np.abs(zeros)


def fibonacci_chain_positions(n_iter=25, tau=None):
    """
    Cut-and-project Fibonacci quasicrystal atom positions via substitution
    L -> LS, S -> L, with length(L) = tau (golden ratio), length(S) = 1.
    Produces a DETERMINISTIC sequence with exactly two distinct gap sizes
    (L and S in ratio tau:1) — a clean structural negative example: no
    level repulsion, no GUE-like statistics expected by construction.
    """
    if tau is None:
        tau = (1 + math.sqrt(5)) / 2
    word = "L"
    for _ in range(n_iter):
        word = word.replace("L", "LX").replace("S", "L").replace("X", "S")
        if len(word) > 2_000_000:
            break
    lengths = np.where(np.array(list(word)) == "L", tau, 1.0)
    positions = np.concatenate([[0.0], np.cumsum(lengths)])
    return positions


def main():
    candidates = {}
    print("Circular billiard spectrum (Bessel zeros^2, integrable -> expect Poisson)...")
    candidates["circular_billiard_bessel2_m200k200"] = circular_billiard_spectrum(200, 200)

    print("Airy function zeros...")
    candidates["airy_zeros_n5000"] = airy_zeros(5000)

    print("Fibonacci quasicrystal atom positions...")
    candidates["fibonacci_chain_positions"] = fibonacci_chain_positions(25)

    results = {}
    for name, seq in candidates.items():
        seq = np.asarray(seq, dtype=float)
        seq = seq[np.isfinite(seq)]
        r = lib.mean_ratio(seq)
        results[name] = {"n_raw": len(seq), "mean_ratio": r}
        print(lib.summarize({"name": name, "n": len(seq), "mean_ratio": r}))

    outpath = os.path.join(RESULTS_DIR, "group_bc_local_ratios.json")
    with open(outpath, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved to {outpath}")

    raw_out = {name: candidates[name].tolist() for name in candidates}
    with open(os.path.join(RESULTS_DIR, "group_bc_local_raw.json"), "w") as f:
        json.dump(raw_out, f)

    return results


if __name__ == "__main__":
    main()
