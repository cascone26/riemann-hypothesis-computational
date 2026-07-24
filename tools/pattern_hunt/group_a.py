"""
Group A candidates — cheap, local, likely-unchecked-against-GUE sequences.

Run: python3 group_a.py
Produces results/pattern_hunt/group_a.json + printed summary table.
"""
import os
import sys
import json
import math
import numpy as np
import networkx as nx

sys.path.insert(0, os.path.dirname(__file__))
import lib

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
                            "results", "pattern_hunt")
os.makedirs(RESULTS_DIR, exist_ok=True)


# =============================================================================
# Candidate generators
# =============================================================================

def prime_gaps(limit=2_000_000):
    """Known NON-GUE (Poisson-ish, actually sub-Poisson at short range due to
    Hardy-Littlewood correlations) — good negative sanity example."""
    from sympy import sieve
    primes = list(sieve.primerange(2, limit))
    return np.array(primes, dtype=float)


def graph_spectrum(G, kind="laplacian"):
    if kind == "laplacian":
        M = nx.laplacian_matrix(G).toarray().astype(float)
    else:
        M = nx.adjacency_matrix(G).toarray().astype(float)
    eigs = np.linalg.eigvalsh(M)
    return eigs


def orthopoly_roots():
    n = 3000
    out = {}
    herm_x, _ = np.polynomial.hermite.hermgauss(n)
    out["hermite_roots_n3000"] = herm_x
    leg_x, _ = np.polynomial.legendre.leggauss(n)
    out["legendre_roots_n3000"] = leg_x
    lag_x, _ = np.polynomial.laguerre.laggauss(n)
    out["laguerre_roots_n3000"] = lag_x
    k = np.arange(1, n + 1)
    cheb_x = np.cos((2 * k - 1) * math.pi / (2 * n))
    out["chebyshev_roots_n3000"] = cheb_x
    return out


def sum_of_two_squares(limit=500_000):
    """Ideal norms of Z[i] (Gaussian integers) — n represented iff every
    prime = 3 (mod 4) in its factorization occurs to an even power."""
    is_valid = np.ones(limit + 1, dtype=bool)
    is_valid[0] = False
    sieve_factor = np.arange(limit + 1)
    n = 2
    temp = np.arange(limit + 1)
    # factor via smallest prime factor sieve
    spf = np.zeros(limit + 1, dtype=np.int64)
    for i in range(2, limit + 1):
        if spf[i] == 0:
            for j in range(i, limit + 1, i):
                if spf[j] == 0:
                    spf[j] = i
    good = []
    for m in range(1, limit + 1):
        x = m
        ok = True
        while x > 1:
            p = spf[x]
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            if p % 4 == 3 and e % 2 == 1:
                ok = False
                break
        if ok:
            good.append(m)
    return np.array(good, dtype=float)


def ideal_norms_quadratic_form(a, b, c, bound=200000, max_count=200000):
    """Norm-form values a*x^2 + b*x*y + c*y^2 for a fixed binary quadratic
    form (ideal norms of the corresponding quadratic order), sorted, distinct."""
    vals = set()
    lim = int(math.isqrt(bound // max(1, min(a, c)))) + 2
    for x in range(-lim, lim + 1):
        for y in range(-lim, lim + 1):
            if x == 0 and y == 0:
                continue
            v = a * x * x + b * x * y + c * y * y
            if 0 < v <= bound:
                vals.add(v)
    arr = np.array(sorted(vals), dtype=float)
    return arr[:max_count]


def cf_convergent_denominators(alpha_name, alpha, n_terms=2000):
    """log(convergent denominator) of a real number's continued fraction —
    denominators grow geometrically (Levy's constant), so we test the
    log-sequence (arbitrary-precision recurrence via Python bigints) rather
    than the raw denominators which overflow float64 after ~700 terms.
    Uses mpmath at high precision (dps scaled to n_terms) so the CF
    expansion itself stays valid for all n_terms requested — a naive
    float64 CF iteration only gives ~15-20 valid terms before rounding
    error takes over. Deterministic/periodic for quadratic irrationals
    (sqrt2, golden) — expect NOT GUE (negative example). Non-periodic for
    pi/e — genuinely unclear a priori, worth the cheap check."""
    import mpmath
    mpmath.mp.dps = max(50, n_terms * 2)
    x = mpmath.mpf(alpha) if not isinstance(alpha, str) else {
        "pi": mpmath.pi, "e": mpmath.e, "sqrt2": mpmath.sqrt(2),
        "golden": (1 + mpmath.sqrt(5)) / 2,
    }[alpha]
    a0 = int(mpmath.floor(x))
    coeffs = [a0]
    frac = x - a0
    for _ in range(n_terms - 1):
        if frac < mpmath.mpf(10) ** (-(mpmath.mp.dps - 10)):
            break
        frac = 1 / frac
        ai = int(mpmath.floor(frac))
        coeffs.append(ai)
        frac = frac - ai
    h_prev2, h_prev1 = 1, coeffs[0]
    k_prev2, k_prev1 = 0, 1
    denoms = [1, k_prev1]
    for a in coeffs[1:]:
        h = a * h_prev1 + h_prev2
        k = a * k_prev1 + k_prev2
        denoms.append(k)
        h_prev2, h_prev1 = h_prev1, h
        k_prev2, k_prev1 = k_prev1, k
    log_denoms = [math.log(d) for d in denoms[1:] if d > 0]
    return np.array(log_denoms, dtype=float)


# =============================================================================
# Runner
# =============================================================================

def main():
    candidates = {}

    print("Building candidates...")

    print("  prime gaps (positions, up to 2M)...")
    candidates["prime_gaps_positions_2M"] = prime_gaps(2_000_000)

    print("  random regular graph spectra...")
    for d in [3, 4, 6, 10]:
        for N in [1200]:
            G = nx.random_regular_graph(d, N, seed=1)
            for kind in ["laplacian", "adjacency"]:
                eigs = graph_spectrum(G, kind)
                candidates[f"regular_graph_d{d}_N{N}_{kind}"] = eigs

    print("  Erdos-Renyi graph spectra...")
    for p in [0.01, 0.05, 0.2]:
        G = nx.erdos_renyi_graph(1200, p, seed=2)
        for kind in ["laplacian", "adjacency"]:
            eigs = graph_spectrum(G, kind)
            candidates[f"erdos_renyi_p{p}_N1200_{kind}"] = eigs

    print("  Watts-Strogatz small-world spectra...")
    G = nx.watts_strogatz_graph(1200, 6, 0.3, seed=3)
    for kind in ["laplacian", "adjacency"]:
        candidates[f"watts_strogatz_N1200_{kind}"] = graph_spectrum(G, kind)

    print("  Barabasi-Albert scale-free spectra...")
    G = nx.barabasi_albert_graph(1200, 4, seed=4)
    for kind in ["laplacian", "adjacency"]:
        candidates[f"barabasi_albert_N1200_{kind}"] = graph_spectrum(G, kind)

    print("  random bipartite regular graph...")
    G = nx.bipartite.random_graph(600, 600, 0.02, seed=5)
    for kind in ["laplacian", "adjacency"]:
        candidates[f"bipartite_random_N1200_{kind}"] = graph_spectrum(G, kind)

    print("  orthogonal polynomial roots...")
    candidates.update(orthopoly_roots())

    print("  sum-of-two-squares (Z[i] ideal norms)...")
    candidates["ideal_norms_Zi_500k"] = sum_of_two_squares(500_000)

    print("  ideal norms of other imaginary quadratic orders...")
    # discriminant -3 (Eisenstein): x^2+xy+y^2 ; disc -7: x^2+xy+2y^2 ; disc -8: x^2+2y^2
    candidates["ideal_norms_disc-3_Eisenstein"] = ideal_norms_quadratic_form(1, 1, 1, bound=500000)
    candidates["ideal_norms_disc-7"] = ideal_norms_quadratic_form(1, 1, 2, bound=500000)
    candidates["ideal_norms_disc-8"] = ideal_norms_quadratic_form(1, 0, 2, bound=500000)
    candidates["ideal_norms_disc-20"] = ideal_norms_quadratic_form(1, 0, 5, bound=500000)

    print("  continued fraction convergent denominators...")
    candidates["cf_log_denoms_sqrt2"] = cf_convergent_denominators("sqrt2", "sqrt2")
    candidates["cf_log_denoms_golden"] = cf_convergent_denominators("golden", "golden")
    candidates["cf_log_denoms_pi"] = cf_convergent_denominators("pi", "pi", n_terms=3000)
    candidates["cf_log_denoms_e"] = cf_convergent_denominators("e", "e", n_terms=3000)

    print(f"\n{len(candidates)} candidates built. Running <r> sweep...\n")

    results = {}
    for name, seq in candidates.items():
        seq = np.asarray(seq, dtype=float)
        seq = seq[np.isfinite(seq)]
        r = lib.mean_ratio(seq)
        results[name] = {"n_raw": len(seq), "mean_ratio": r}
        flag = " <-- GUE-like!" if r["closest"] == "GUE" and abs(r.get("sigma_from_gue", 99)) < 5 else ""
        print(f"  {name:42s} N={len(seq):>7d}  <r>={r['mean']:.4f}+/-{r['stderr']:.4f}  closest={r['closest']:8s}{flag}")

    outpath = os.path.join(RESULTS_DIR, "group_a_ratios.json")
    with open(outpath, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved to {outpath}")

    # Save raw sequences for candidates that look interesting (closest==GUE with n>=200)
    interesting = {name: candidates[name].tolist() for name, r in
                    [(n, results[n]["mean_ratio"]) for n in results]
                    if r["closest"] == "GUE" and r["n"] >= 200}
    # (the comprehension above is written oddly to reuse r var name; rebuild cleanly)
    interesting = {}
    for name, res in results.items():
        r = res["mean_ratio"]
        if r["closest"] == "GUE" and r["n"] >= 200:
            interesting[name] = candidates[name].tolist()
    with open(os.path.join(RESULTS_DIR, "group_a_gue_hits_raw.json"), "w") as f:
        json.dump(interesting, f)
    print(f"{len(interesting)} candidates flagged closest=GUE with N>=200 -> group_a_gue_hits_raw.json")

    return results


if __name__ == "__main__":
    main()
