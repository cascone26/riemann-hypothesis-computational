"""
Connes Consensus: Multi-Lambda Zero Filtering

Spurious zeros in ξ̂(z) are artifacts of the Fourier grid at z = 2πj/L.
Since L = log(λ²), different lambda values put the grid poles at different
positions. Real zeta zeros are constant; spurious zeros shift with λ.

Strategy:
  - Run the Connes spectral triple for λ² ∈ {14, 30, 50}
  - Find zeros of ξ̂ for each
  - Take consensus: only keep zeros that appear in ≥2 of the 3 runs
    (within a tolerance window)

Expected outcome: spurious zeros eliminated, hit rate improves significantly.
"""

import mpmath
import numpy as np
from scipy.optimize import brentq
import json
import os
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR    = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
PLOTS_DIR   = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")

mpmath.mp.dps = 80
N_TRUNC = 100   # slightly smaller than 120 for speed × 3 runs


def load_zeros(count=100):
    path = os.path.join(DATA_DIR, "zeros_100k.txt")
    zeros = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i >= count:
                break
            line = line.strip()
            if line:
                zeros.append(float(line))
    return np.array(zeros)


def primes_up_to(n):
    sieve = list(range(n+1))
    for i in range(2, int(n**0.5)+1):
        if sieve[i]:
            for j in range(i*i, n+1, i):
                sieve[j] = 0
    return [x for x in sieve[2:] if x]


def von_mangoldt(k, primes):
    if k <= 1:
        return 0.0
    for p in primes:
        if p > k:
            break
        pk = p
        while pk <= k:
            if pk == k:
                return float(mpmath.log(p))
            pk *= p
    return 0.0


def W_archimedean(n, m, L):
    pi = float(mpmath.pi)
    if n == m:
        if n == 0:
            return float(-mpmath.euler - mpmath.log(2))
        s = pi * n / L
        val = float(mpmath.re(mpmath.digamma(mpmath.mpf('0.25') + 1j*s)))
        return val + float(mpmath.log(2))
    else:
        k = n - m
        s_avg = pi * (abs(n) + abs(m)) / (2*L)
        if s_avg < 0.01:
            s_avg = 0.01
        psi_im = float(mpmath.im(mpmath.digamma(mpmath.mpf('0.25') + 1j*s_avg)))
        return ((-1)**abs(k)) * psi_im / (pi * abs(k))


def W_boundary(n, m, L):
    pi = float(mpmath.pi)
    L_mp = mpmath.mpf(str(L))
    num = 32*L_mp * mpmath.sinh(L_mp/4)**2 * (L_mp**2 - 16*pi**2*m*n)
    den = (L_mp**2 + 16*pi**2*m**2) * (L_mp**2 + 16*pi**2*n**2)
    if abs(float(den)) < 1e-50:
        return 0.0
    return float(num / den)


def W_prime(n, m, L, primes):
    pi = float(mpmath.pi)
    k_max = int(float(mpmath.exp(mpmath.mpf(str(L))))) + 1
    total = 0.0
    for k in range(2, k_max + 1):
        lam_k = von_mangoldt(k, primes)
        if lam_k == 0:
            continue
        y = float(mpmath.log(k))
        weight = lam_k * k**(-0.5)
        if n != m:
            q = (np.sin(2*pi*m*y/L) - np.sin(2*pi*n*y/L)) / (pi * (n - m))
        else:
            q = 2*(1 - abs(y)/L) * np.cos(2*pi*n*y/L)
        total += weight * q
    return total


def run_connes(lambda_sq, N=N_TRUNC, verbose=True):
    """Full Connes run for a given lambda²."""
    lam = float(mpmath.sqrt(lambda_sq))
    L   = float(2 * mpmath.log(mpmath.sqrt(lambda_sq)))
    primes = primes_up_to(lambda_sq)
    dim = 2*N + 1

    if verbose:
        print(f"\n  λ²={lambda_sq}: L={L:.4f}, primes={primes}, matrix {dim}×{dim}")

    # Build matrix
    tau = np.zeros((dim, dim))
    t0 = time.time()
    for i in range(dim):
        n = i - N
        for j in range(i, dim):
            m = j - N
            val = (W_archimedean(n, m, L)
                   + W_boundary(n, m, L)
                   - W_prime(n, m, L, primes))
            tau[i,j] = val
            tau[j,i] = val
    if verbose:
        print(f"  Matrix built in {time.time()-t0:.1f}s")

    # Min eigenvector
    evals, evecs = np.linalg.eigh(tau)
    eps = evals[0]
    xi  = evecs[:, 0].copy()
    delta = xi.sum() / np.sqrt(L)
    if abs(delta) > 1e-12:
        xi /= delta

    if verbose:
        print(f"  ε_N = {eps:.6f}")

    # Find zeros of ξ̂
    poles = np.array([2.0*np.pi*j/L for j in range(-N, N+1)])

    def xi_hat(z):
        prefix = (2.0/np.sqrt(L)) * np.sin(z*L/2.0)
        diffs  = z - poles
        mask   = np.abs(diffs) > 1e-12
        return prefix * np.sum(xi[mask] / diffs[mask])

    z_grid = np.linspace(12.0, 200.0, 60000)
    vals   = np.array([xi_hat(z) for z in z_grid])

    found = []
    for i in range(len(z_grid)-1):
        if vals[i]*vals[i+1] < 0 and not (np.isnan(vals[i]) or np.isnan(vals[i+1])):
            try:
                root = brentq(xi_hat, z_grid[i], z_grid[i+1], xtol=1e-8)
                found.append(root)
            except Exception:
                pass

    if verbose:
        print(f"  Found {len(found)} zeros (z > 12)")

    return {
        "lambda_sq": lambda_sq, "L": L, "primes": primes,
        "epsilon_N": float(eps), "zeros": found,
    }


def consensus(runs, tol=0.8):
    """
    Keep only zeros that appear in ≥2 of the runs within `tol`.
    For each consensus zero, take the median position.
    """
    all_zeros = []
    for r in runs:
        all_zeros.extend([(z, r["lambda_sq"]) for z in r["zeros"]])
    all_zeros.sort(key=lambda x: x[0])

    used = [False] * len(all_zeros)
    groups = []

    for i, (z0, lsq0) in enumerate(all_zeros):
        if used[i]:
            continue
        group = [(z0, lsq0)]
        for j, (z1, lsq1) in enumerate(all_zeros[i+1:], start=i+1):
            if z1 - z0 > tol:
                break
            if not used[j] and lsq1 != lsq0:
                group.append((z1, lsq1))
                used[j] = True
        used[i] = True
        # Need ≥2 different lambdas in the group
        lambdas = set(lsq for _, lsq in group)
        if len(lambdas) >= 2:
            groups.append(group)

    consensus_zeros = [np.median([z for z, _ in g]) for g in groups]
    return sorted(consensus_zeros), groups


def main():
    print("CONNES CONSENSUS: Multi-Lambda Zero Filtering")
    print("=" * 55)
    print(f"N={N_TRUNC}, λ² values: 14, 30, 50")
    print()

    ref_zeros = load_zeros(100)
    lambda_sqs = [14, 30, 50]
    runs = []

    for lsq in lambda_sqs:
        r = run_connes(lsq, N=N_TRUNC, verbose=True)
        runs.append(r)

    print("\n" + "=" * 55)
    print("CONSENSUS FILTERING")
    print("=" * 55)

    con_zeros, groups = consensus(runs, tol=1.0)
    print(f"Consensus zeros (appear in ≥2 lambdas): {len(con_zeros)}")
    print()

    # Match against reference
    matched_tight = 0
    matched_loose = 0
    comparisons = []
    print(f"{'#':>4}  {'Consensus':>10}  {'Nearest ζ':>10}  {'Error':>10}  {'Match?':>7}")
    print("-" * 52)
    for i, cz in enumerate(con_zeros[:60]):
        diffs = np.abs(ref_zeros - cz)
        idx   = np.argmin(diffs)
        err   = diffs[idx]
        tag   = "✓" if err < 0.5 else ("~" if err < 2.0 else " ")
        if err < 0.5: matched_tight += 1
        if err < 2.0: matched_loose += 1
        comparisons.append({"consensus": cz, "zero": ref_zeros[idx], "error": float(err)})
        print(f"{i+1:>4}  {cz:>10.4f}  {ref_zeros[idx]:>10.4f}  {err:>10.6f}  {tag:>7}")

    n_shown = min(len(con_zeros), 60)
    print()
    print(f"Summary ({n_shown} consensus zeros):")
    print(f"  Matched tight (<0.5): {matched_tight}/{n_shown}")
    print(f"  Matched loose (<2.0): {matched_loose}/{n_shown}")
    print()
    print(f"v1  N=30,  1 lambda:    12/20 (60%)")
    print(f"v2  N=120, 1 lambda:    28/55 (51%, but 3× more zeros covered)")
    print(f"v3  N={N_TRUNC}, 3 lambdas: {matched_tight}/{n_shown} ({100*matched_tight/n_shown:.0f}%)")

    # Per-lambda zero counts for comparison
    print()
    print("Per-lambda stats:")
    for r in runs:
        lz = np.array(r["zeros"])
        tight = sum(1 for z in lz if any(abs(z-rz) < 0.5 for rz in ref_zeros[:50]))
        print(f"  λ²={r['lambda_sq']}: {len(lz)} found, {tight} match tight")

    # ── Plots ────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 1, figsize=(16, 8))

    ax = axes[0]
    # Show all per-lambda zeros as ticks, consensus as filled circles
    colors_lsq = {14: 'blue', 30: 'orange', 50: 'green'}
    for r in runs:
        for z in r["zeros"]:
            ax.axvline(z, color=colors_lsq[r["lambda_sq"]], alpha=0.15, linewidth=0.6)
    for z in ref_zeros[:50]:
        ax.axvline(z, color='red', alpha=0.5, linewidth=1.0, linestyle=':')
    ax.scatter(con_zeros[:60], [0.5]*len(con_zeros[:60]), c='black',
               zorder=5, s=20, label='Consensus zeros')
    ax.set_xlim(12, 150)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel('z')
    ax.set_title('Connes Consensus: per-λ zeros (blue/orange/green), ζ zeros (red), consensus (black)')
    ax.legend(handles=[
        plt.Line2D([0],[0], color='blue',   label='λ²=14'),
        plt.Line2D([0],[0], color='orange', label='λ²=30'),
        plt.Line2D([0],[0], color='green',  label='λ²=50'),
        plt.Line2D([0],[0], color='red',    linestyle=':', label='ζ zeros'),
        plt.scatter([],[], c='black', label='Consensus'),
    ], fontsize=7)

    ax = axes[1]
    if comparisons:
        errs   = [c["error"] for c in comparisons]
        idxs   = list(range(1, len(errs)+1))
        colors = ['green' if e < 0.5 else ('orange' if e < 2.0 else 'red') for e in errs]
        ax.bar(idxs, errs, color=colors, alpha=0.8)
        ax.axhline(0.5, color='green',  linestyle='--', linewidth=1, label='tight (0.5)')
        ax.axhline(2.0, color='orange', linestyle='--', linewidth=1, label='loose (2.0)')
        ax.set_xlabel('Consensus zero index')
        ax.set_ylabel('|consensus − nearest ζ zero|')
        ax.set_title(f'Error per consensus zero — {matched_tight}/{n_shown} tight, {matched_loose}/{n_shown} loose')
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.suptitle(
        f'Connes Consensus (N={N_TRUNC}, λ²∈{{14,30,50}}) — '
        f'{matched_tight}/{n_shown} tight matches ({100*matched_tight/n_shown:.0f}%)',
        fontsize=11
    )
    plt.tight_layout()
    out = os.path.join(PLOTS_DIR, "connes_consensus.png")
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"\nPlot saved: {out}")

    results = {
        "N_TRUNC": N_TRUNC,
        "lambda_sqs": lambda_sqs,
        "consensus_zeros": con_zeros[:60],
        "matched_tight": matched_tight,
        "matched_loose": matched_loose,
        "n_shown": n_shown,
        "per_lambda": [{"lambda_sq": r["lambda_sq"], "n_found": len(r["zeros"]),
                        "epsilon_N": r["epsilon_N"]} for r in runs],
        "comparisons": comparisons[:60],
    }
    with open(os.path.join(RESULTS_DIR, "connes_consensus.json"), "w") as f:
        json.dump(results, f, indent=2)
    print("Results saved.")


if __name__ == "__main__":
    main()
