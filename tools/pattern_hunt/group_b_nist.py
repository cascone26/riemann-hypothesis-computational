"""
Group B — NIST Atomic Spectra Database energy levels (heavy/complex atoms
with dense level structure: Th I, Th II, Nd II). Classic quantum-chaos test
case (complex atoms are known from decades of atomic-physics literature to
show GOE-like level statistics WITHIN a fixed J/parity block — mixing
different J blocks is a well-known pitfall that spuriously looks Poisson
since crossings between different-symmetry sequences aren't repelled).

We report BOTH: naive mixed-J statistics (expected artificially Poisson-ish)
and same-J-block statistics (the methodologically correct test) for honesty.
"""
import os
import sys
import re
import json
import numpy as np
from collections import defaultdict

sys.path.insert(0, os.path.dirname(__file__))
import lib

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
                            "results", "pattern_hunt")


def parse_nist_csv(path):
    levels_by_j = defaultdict(list)
    all_levels = []
    with open(path, encoding="latin-1") as f:
        next(f)  # header
        for line in f:
            fields = re.findall(r'"=""([^"]*)"""', line)
            if len(fields) < 5:
                continue
            conf, term, j, prefix, level = fields[0], fields[1], fields[2], fields[3], fields[4]
            try:
                lvl = float(level)
            except ValueError:
                continue
            if lvl <= 0 and len(all_levels) > 0:
                continue  # skip duplicate ground states etc
            all_levels.append(lvl)
            levels_by_j[j].append(lvl)
    return all_levels, levels_by_j


def main():
    files = {
        "Th_I": "/tmp/nist_Th_I.csv",
        "Th_II": "/tmp/nist_Th_II.csv",
        "Nd_II": "/tmp/nist_Nd_II.csv",
    }
    results = {}
    for name, path in files.items():
        if not os.path.exists(path):
            continue
        all_levels, by_j = parse_nist_csv(path)
        all_levels = sorted(set(all_levels))
        r_mixed = lib.mean_ratio(all_levels)
        print(f"{name} MIXED-J  N={len(all_levels)}  <r>={r_mixed['mean']:.4f}+/-{r_mixed['stderr']:.4f}  closest={r_mixed['closest']}")
        results[f"{name}_mixed_J"] = {"n": len(all_levels), "mean_ratio": r_mixed}

        # best single-J block
        best_j, best_levels = max(by_j.items(), key=lambda kv: len(set(kv[1])))
        best_levels = sorted(set(best_levels))
        r_j = lib.mean_ratio(best_levels)
        print(f"{name} J={best_j:>5s} (largest block) N={len(best_levels)}  <r>={r_j['mean']:.4f}+/-{r_j['stderr']:.4f}  closest={r_j['closest']}")
        results[f"{name}_J{best_j}_block"] = {"n": len(best_levels), "j": best_j, "mean_ratio": r_j}

        # show the J-block size distribution for transparency
        sizes = sorted(((j, len(set(v))) for j, v in by_j.items()), key=lambda x: -x[1])
        print(f"   J-block sizes: {sizes[:6]}")

    outpath = os.path.join(RESULTS_DIR, "group_b_nist_ratios.json")
    with open(outpath, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved to {outpath}")
    return results


if __name__ == "__main__":
    main()
