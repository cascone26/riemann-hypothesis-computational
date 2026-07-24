"""Positive control (zeta zeros) and negative control (Poisson/random) —
run BEFORE trusting any new candidate result."""

import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import lib

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def load_zeta(n=20000):
    path = os.path.join(DATA_DIR, "zeros_odlyzko_100k")
    zeros = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i >= n:
                break
            line = line.strip()
            if line:
                zeros.append(float(line))
    return zeros


def zeta_unfold(zeros):
    """Zeta zeros have a KNOWN smooth trend (Riemann-von Mangoldt); use it
    directly rather than generic unfolding, matching exact_gue.py."""
    import math
    t = np.asarray(zeros)
    twopi = 2.0 * np.pi
    N = (t / twopi) * np.log(t / twopi) - t / twopi + 7.0 / 8.0
    return np.diff(N)


def main():
    print("=" * 70)
    print("CONTROLS")
    print("=" * 70)

    zeros = load_zeta(100000)
    r_zeta = lib.mean_ratio(zeros)
    print(f"\nPOSITIVE CONTROL — zeta zeros (N={len(zeros)}):")
    print(f"  <r> = {r_zeta['mean']:.5f} +/- {r_zeta['stderr']:.5f}  (expect ~0.6027, GUE)")

    unfolded = zeta_unfold(zeros)
    centers, pdf = lib.spacing_pdf(unfolded)
    ks = lib.ks_against_ensembles(unfolded)
    pc_centers, pc_emp, pc_gue, N = lib.pair_correlation(unfolded)
    pc_mad = float(np.mean(np.abs(pc_emp[5:] - pc_gue[5:])))
    nv = lib.number_variance(unfolded)
    print(f"  unfolded mean={np.mean(unfolded):.4f} var={np.var(unfolded):.4f} (GUE var ~0.178)")
    print(f"  KS D-stats: " + ", ".join(f"{k}={v['D']:.4f}" for k, v in ks.items()))
    print(f"  pair-corr MAD vs GUE (x>0.125): {pc_mad:.4f}")
    for row in nv:
        print(f"  Sigma^2(L={row['L']:>5}): emp={row['empirical_variance']:.4f} gue={row['gue_variance']:.4f} poisson={row['poisson_variance']:.4f}")

    np.random.seed(42)
    random_seq = np.sort(np.random.uniform(0, 20000, 20000))
    r_rand = lib.mean_ratio(random_seq)
    print(f"\nNEGATIVE CONTROL — uniform random N={len(random_seq)}):")
    print(f"  <r> = {r_rand['mean']:.5f} +/- {r_rand['stderr']:.5f}  (expect ~0.3863, Poisson)")

    ru = lib.unfold_local(random_seq)
    ks_r = lib.ks_against_ensembles(ru)
    print(f"  KS D-stats: " + ", ".join(f"{k}={v['D']:.4f}" for k, v in ks_r.items()))

    # Note: convergence of finite-height zeta <r> to the asymptotic GUE value
    # 0.6027 is logarithmically slow (documented in STATUS.md — 100K zeros
    # gives 0.6109, matching nn_ratio_exact.py's prior measurement exactly).
    # Validate on ensemble ranking + reasonable tolerance, not tight match.
    ok = (r_zeta["closest"] == "GUE" and abs(r_zeta["mean"] - lib.GUE_R) < 0.02
          and r_rand["closest"] == "Poisson" and abs(r_rand["mean"] - lib.POISSON_R) < 0.02)
    print(f"\nPIPELINE VALIDATED: {ok}")

    out = {
        "zeta": {"mean_ratio": r_zeta, "unfolded_var": float(np.var(unfolded)),
                 "ks": ks, "pair_corr_mad": pc_mad, "number_variance": nv},
        "random": {"mean_ratio": r_rand, "ks": ks_r},
        "pipeline_validated": bool(ok),
    }
    outdir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
                           "results", "pattern_hunt")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "controls.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved to {os.path.join(outdir, 'controls.json')}")
    return out


if __name__ == "__main__":
    main()
