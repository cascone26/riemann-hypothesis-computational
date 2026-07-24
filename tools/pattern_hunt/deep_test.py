"""
Deep-test a named candidate sequence: full suite (P(s), pair correlation,
Delta_3/number variance) + plot. Only run this on candidates whose <r> was
already close to the GUE value (0.6027) in the cheap first pass.
"""
import os
import sys
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
import lib

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
                            "results", "pattern_hunt")
PLOTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
                          "plots", "pattern_hunt")
os.makedirs(PLOTS_DIR, exist_ok=True)


def deep_test(name, seq, unfold_method="local", notes=""):
    result = lib.run_full_suite(name, seq, unfold_method=unfold_method, notes=notes)
    print(lib.summarize(result))
    if "spacing_pdf" not in result:
        return result

    centers = np.array(result["spacing_pdf"]["centers"])
    pdf = np.array(result["spacing_pdf"]["pdf"])

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

    ax = axes[0]
    ax.bar(centers, pdf, width=centers[1] - centers[0], alpha=0.5, color="steelblue", label=name)
    s_grid = np.linspace(0.001, centers[-1], 300)
    ax.plot(s_grid, lib.gue_wigner_pdf(s_grid), "r-", lw=2, label="GUE")
    ax.plot(s_grid, lib.goe_wigner_pdf(s_grid), "g--", lw=1.5, label="GOE")
    ax.plot(s_grid, lib.poisson_pdf(s_grid), "k:", lw=1.5, label="Poisson")
    ax.set_title(f"P(s) — {name}")
    ax.set_xlabel("s")
    ax.legend(fontsize=8)

    ax = axes[1]
    if unfold_method == "poly":
        unfolded = lib.unfold_poly(seq)
    elif unfold_method == "rank":
        unfolded = lib.unfold_rank(seq)
    else:
        unfolded = lib.unfold_local(seq)
    pc_centers, pc_emp, pc_gue, N = lib.pair_correlation(unfolded)
    ax.plot(pc_centers, pc_emp, "b-", lw=0.8, label="empirical")
    ax.plot(pc_centers, pc_gue, "r--", lw=1.5, label="GUE 1-(sinc)^2")
    ax.set_title(f"Pair correlation R2(x)")
    ax.set_xlabel("x")
    ax.legend(fontsize=8)

    ax = axes[2]
    nv = result["number_variance"]
    Ls = [row["L"] for row in nv]
    emp = [row["empirical_variance"] for row in nv]
    gue = [row["gue_variance"] for row in nv]
    poi = [row["poisson_variance"] for row in nv]
    ax.plot(Ls, emp, "bo-", label="empirical")
    ax.plot(Ls, gue, "r--", label="GUE")
    ax.plot(Ls, poi, "k:", label="Poisson")
    ax.set_title("Number variance Sigma^2(L)")
    ax.set_xlabel("L")
    ax.legend(fontsize=8)

    plt.tight_layout()
    safe = "".join(c if c.isalnum() else "_" for c in name)
    outpath = os.path.join(PLOTS_DIR, f"{safe}.png")
    plt.savefig(outpath, dpi=130)
    plt.close(fig)
    result["plot"] = outpath
    print(f"  plot -> {outpath}")
    return result


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--from-json", help="load {name: [values...]} dict and deep-test all")
    parser.add_argument("--unfold", default="local", choices=["local", "poly", "rank"])
    args = parser.parse_args()

    all_results = {}
    if args.from_json:
        with open(args.from_json) as f:
            data = json.load(f)
        for name, seq in data.items():
            all_results[name] = deep_test(name, np.array(seq), unfold_method=args.unfold)

    outpath = os.path.join(RESULTS_DIR, "deep_test_results.json")
    with open(outpath, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nSaved to {outpath}")
