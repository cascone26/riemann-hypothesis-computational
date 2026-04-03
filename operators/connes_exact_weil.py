"""
Connes Exact Weil Matrix from Known Zeros

Instead of approximating the archimedean off-diagonal terms W_R(V_n, V_m),
we compute the Weil matrix directly from known zeros:

    τ(V_n, V_m) = Σ_{γ: ζ(1/2+iγ)=0} ê_n(γ) · ê_m(γ)*

where ê_n(γ) is the Mellin transform of V_n evaluated on the critical line:

    ê_n(γ) = (1/√L) · e^{i(γ + 2πn/L)·L/2} · 2·sin((γ + 2πn/L)·L/2) / (γ + 2πn/L)

This is the exact Weil explicit formula — no approximation of the archimedean term.
We add the prime and boundary contributions exactly as before.

The minimum eigenvector of this exact matrix should have zeros that ARE the
remaining zeta zeros (those not in our input set), if RH is true.

Test: use first K zeros to build the matrix, find minimum eigenvector,
check if its zeros match zeros K+1, K+2, ...
"""

import numpy as np
import mpmath
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

N_TRUNC   = 80      # matrix size: 161×161
LAMBDA_SQ = 14
L         = float(2 * mpmath.log(mpmath.sqrt(LAMBDA_SQ)))
PRIMES    = [2, 3, 5, 7, 11, 13]
K_ZEROS   = 500     # how many zeros to use for building the Weil matrix


def load_zeros(count=2000):
    path = os.path.join(DATA_DIR, "zeros_odlyzko_2M")
    try:
        zeros = []
        with open(path) as f:
            for i, line in enumerate(f):
                if i >= count:
                    break
                line = line.strip()
                if line:
                    zeros.append(float(line))
        if len(zeros) >= count:
            return np.array(zeros)
    except Exception:
        pass
    # fallback to 100k file
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


def e_hat_n(gamma, n, L):
    """
    Mellin transform of V_n evaluated at γ:
    ê_n(γ) = (1/√L) · e^{i(γ+2πn/L)L/2} · 2sin((γ+2πn/L)L/2) / (γ+2πn/L)

    Returns complex float.
    """
    freq = gamma + 2 * np.pi * n / L
    if abs(freq) < 1e-12:
        # L'Hopital: limit is (1/√L) · e^0 · L = √L
        return np.sqrt(L) + 0j
    half_arg = freq * L / 2.0
    return (1.0 / np.sqrt(L)) * np.exp(1j * half_arg) * 2 * np.sin(half_arg) / freq


def build_weil_matrix_from_zeros(zeros_subset, verbose=True):
    """
    Build τ = Σ_γ ê_n(γ)·ê_m(γ)* + boundary - primes.

    The zero-sum term is the exact archimedean contribution.
    """
    dim = 2 * N_TRUNC + 1
    K   = len(zeros_subset)
    t0  = time.time()

    if verbose:
        print(f"Building {dim}×{dim} exact Weil matrix from K={K} zeros...")

    # ── Zero-sum term (archimedean) ──────────────────────────────
    # Compute ê_n(γ) for all n and all γ at once
    ns     = np.arange(-N_TRUNC, N_TRUNC + 1)   # shape (dim,)
    gammas = zeros_subset                         # shape (K,)

    # e_hat[k, i] = ê_{n_i}(γ_k)   shape: (K, dim)
    freq_matrix  = gammas[:, None] + 2*np.pi*ns[None, :] / L   # (K, dim)
    half_arg     = freq_matrix * L / 2.0
    # Avoid division by zero
    safe_freq    = np.where(np.abs(freq_matrix) < 1e-12, 1.0, freq_matrix)
    e_hat_matrix = ((1/np.sqrt(L))
                    * np.exp(1j * half_arg)
                    * 2 * np.sin(half_arg)
                    / safe_freq)
    # Fix near-zero entries
    near_zero = np.abs(freq_matrix) < 1e-12
    e_hat_matrix[near_zero] = np.sqrt(L)

    # τ_zeros = Σ_γ ê_n · ê_m*  = e_hat^T @ e_hat.conj
    tau_zeros = e_hat_matrix.T @ e_hat_matrix.conj()   # (dim, dim)
    tau = np.real(tau_zeros)   # should be real and symmetric

    if verbose:
        print(f"  Zero-sum term built in {time.time()-t0:.1f}s")

    # ── Boundary term ────────────────────────────────────────────
    pi = np.pi
    for i in range(dim):
        n = i - N_TRUNC
        for j in range(i, dim):
            m = j - N_TRUNC
            num = 32*L * np.sinh(L/4)**2 * (L**2 - 16*pi**2*m*n)
            den = (L**2 + 16*pi**2*m**2) * (L**2 + 16*pi**2*n**2)
            val = num/den if abs(den) > 1e-12 else 0.0
            tau[i, j] += val
            if i != j:
                tau[j, i] += val

    # ── Prime term ───────────────────────────────────────────────
    def von_mangoldt(k):
        for p in PRIMES:
            if p > k: break
            pk = p
            while pk <= k:
                if pk == k: return np.log(p)
                pk *= p
        return 0.0

    k_max = int(np.exp(L)) + 1
    for k in range(2, k_max+1):
        lam_k = von_mangoldt(k)
        if lam_k == 0:
            continue
        y = np.log(k)
        weight = lam_k * k**(-0.5)
        for i in range(dim):
            n = i - N_TRUNC
            for j in range(i, dim):
                m = j - N_TRUNC
                if n != m:
                    q = (np.sin(2*pi*m*y/L) - np.sin(2*pi*n*y/L)) / (pi*(n-m))
                else:
                    q = 2*(1 - abs(y)/L) * np.cos(2*pi*n*y/L)
                tau[i,j] -= weight * q
                if i != j:
                    tau[j,i] -= weight * q

    if verbose:
        print(f"  Full matrix built in {time.time()-t0:.1f}s")

    return tau


def xi_hat_from_xi(xi_vec):
    """Return xi_hat function for the given eigenvector."""
    poles = np.array([2.0*np.pi*j/L for j in range(-N_TRUNC, N_TRUNC+1)])
    def xi_hat(z):
        prefix = (2.0/np.sqrt(L)) * np.sin(z*L/2.0)
        diffs  = z - poles
        mask   = np.abs(diffs) > 1e-12
        return prefix * np.sum(xi_vec[mask] / diffs[mask])
    return xi_hat


def find_zeros(xi_hat_fn, z_lo=12, z_hi=200, n_grid=60000):
    z_grid = np.linspace(z_lo, z_hi, n_grid)
    vals   = np.array([xi_hat_fn(z) for z in z_grid])
    found  = []
    for i in range(len(z_grid)-1):
        if vals[i]*vals[i+1] < 0 and not (np.isnan(vals[i]) or np.isnan(vals[i+1])):
            try:
                root = brentq(xi_hat_fn, z_grid[i], z_grid[i+1], xtol=1e-8)
                found.append(root)
            except Exception:
                pass
    return found


def main():
    print("CONNES EXACT WEIL MATRIX FROM KNOWN ZEROS")
    print("=" * 55)
    print(f"N={N_TRUNC}, λ²={LAMBDA_SQ}, K_input={K_ZEROS} zeros")
    print()

    all_zeros = load_zeros(K_ZEROS + 200)
    print(f"Loaded {len(all_zeros)} zeros. Range: [{all_zeros[0]:.3f}, {all_zeros[-1]:.3f}]")
    print()

    # Split: use first K for matrix, hold out rest for validation
    train_zeros = all_zeros[:K_ZEROS]
    test_zeros  = all_zeros[K_ZEROS:K_ZEROS+200]
    print(f"Training zeros: 1 to {K_ZEROS} (up to γ={train_zeros[-1]:.2f})")
    print(f"Test zeros:     {K_ZEROS+1} to {K_ZEROS+200} (up to γ={test_zeros[-1]:.2f})")
    print()

    # Build exact Weil matrix
    tau = build_weil_matrix_from_zeros(train_zeros, verbose=True)
    print()

    # Minimum eigenvector
    print("Computing eigendecomposition...")
    t0 = time.time()
    evals, evecs = np.linalg.eigh(tau)
    print(f"  Done in {time.time()-t0:.1f}s")
    print(f"  ε_N = {evals[0]:.8f}")
    print(f"  5 smallest: {evals[:5]}")
    print()

    xi = evecs[:, 0].copy()
    delta = xi.sum() / np.sqrt(L)
    if abs(delta) > 1e-12:
        xi /= delta

    xi_hat_fn = xi_hat_from_xi(xi)

    # Find zeros of ξ̂ — looking for zeros BEYOND the training range
    z_lo = train_zeros[-1] + 0.5   # start just after training range
    z_hi = test_zeros[-1] + 5.0

    print(f"Searching for zeros of ξ̂ in [{z_lo:.1f}, {z_hi:.1f}]...")
    t0 = time.time()
    found_beyond = find_zeros(xi_hat_fn, z_lo=z_lo, z_hi=z_hi, n_grid=40000)
    print(f"  Found {len(found_beyond)} zeros beyond training range in {time.time()-t0:.1f}s")
    print()

    # Also find within training range for comparison
    found_within = find_zeros(xi_hat_fn, z_lo=12, z_hi=train_zeros[-1], n_grid=60000)
    print(f"Found {len(found_within)} zeros within training range")
    print()

    # Match beyond-range zeros against test zeros
    print("VALIDATION: Do ξ̂ zeros match test zeros?")
    print(f"{'#':>4}  {'Found':>10}  {'Test zero':>10}  {'Error':>10}  {'Match?':>7}")
    print("-" * 52)
    matched_tight = 0
    matched_loose = 0
    comparisons = []
    for i, fz in enumerate(found_beyond[:60]):
        diffs = np.abs(test_zeros - fz)
        idx   = np.argmin(diffs)
        err   = diffs[idx]
        tag   = "✓" if err < 0.5 else ("~" if err < 2.0 else " ")
        if err < 0.5: matched_tight += 1
        if err < 2.0: matched_loose += 1
        comparisons.append({"found": fz, "zero": test_zeros[idx], "error": float(err)})
        if i < 30:
            print(f"{i+1:>4}  {fz:>10.4f}  {test_zeros[idx]:>10.4f}  {err:>10.6f}  {tag:>7}")

    n_shown = min(len(found_beyond), 60)
    print()
    print(f"Summary (first {n_shown} found beyond training):")
    print(f"  Matched tight (<0.5): {matched_tight}/{n_shown}")
    print(f"  Matched loose (<2.0): {matched_loose}/{n_shown}")
    print()

    # Also check within-range accuracy
    tight_in = sum(1 for fz in found_within
                   if any(abs(fz-z) < 0.5 for z in train_zeros))
    print(f"Within training range: {tight_in}/{len(found_within)} tight matches")
    print()

    print("KEY QUESTION: If ε_N ≈ 0, ξ̂ zeros should be the zeta zeros.")
    print(f"ε_N = {evals[0]:.8f}")
    if evals[0] < 0.1:
        print("ε_N < 0.1: Strong evidence pointing to RH consistency")
    elif evals[0] < 1.0:
        print("ε_N < 1.0: Moderate — more zeros or larger N needed")
    else:
        print("ε_N > 1.0: Matrix truncation artifacts dominate")

    # ── Plot ────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 1, figsize=(16, 8))

    # ξ̂(z) plot over full range
    ax = axes[0]
    z_plot = np.linspace(train_zeros[-1] - 5, min(test_zeros[-1]+5, z_hi), 30000)
    xi_vals = np.clip([xi_hat_fn(z) for z in z_plot], -3, 3)
    ax.plot(z_plot, xi_vals, 'b-', linewidth=0.4)
    ax.axhline(0, color='k', linewidth=0.5)
    for z in test_zeros[:60]:
        ax.axvline(z, color='red', alpha=0.3, linewidth=0.8)
    for fz in found_beyond[:60]:
        ax.axvline(fz, color='green', alpha=0.3, linewidth=0.8, linestyle='--')
    ax.axvline(train_zeros[-1], color='purple', linewidth=2, linestyle='-.',
               label=f'Training cutoff (γ_{K_ZEROS}={train_zeros[-1]:.1f})')
    ax.set_xlim(train_zeros[-1] - 5, min(z_hi, train_zeros[-1] + 60))
    ax.set_ylim(-3, 3)
    ax.set_xlabel('z')
    ax.set_ylabel('ξ̂(z) [clipped ±3]')
    ax.set_title(f'Exact Weil (K={K_ZEROS}) ξ̂ beyond training range — red=ζ zeros, green=predicted')
    ax.legend()
    ax.grid(True, alpha=0.2)

    # Error plot
    ax = axes[1]
    if comparisons:
        errs   = [c["error"] for c in comparisons]
        idxs   = list(range(1, len(errs)+1))
        colors = ['green' if e<0.5 else ('orange' if e<2.0 else 'red') for e in errs]
        ax.bar(idxs, errs, color=colors, alpha=0.8)
        ax.axhline(0.5, color='green',  linestyle='--', linewidth=1, label='tight (0.5)')
        ax.axhline(2.0, color='orange', linestyle='--', linewidth=1, label='loose (2.0)')
        ax.set_xlabel('Found zero index (beyond training)')
        ax.set_ylabel('Error')
        ax.set_title(f'Prediction error — {matched_tight}/{n_shown} tight ({100*matched_tight/max(n_shown,1):.0f}%)')
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.suptitle(
        f'Exact Weil from K={K_ZEROS} zeros, N={N_TRUNC} | ε_N={evals[0]:.6f} | '
        f'Predicted {matched_tight}/{n_shown} tight beyond training',
        fontsize=10
    )
    plt.tight_layout()
    out = os.path.join(PLOTS_DIR, "connes_exact_weil.png")
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"Plot saved: {out}")

    results = {
        "N_TRUNC": N_TRUNC, "K_zeros": K_ZEROS,
        "epsilon_N": float(evals[0]),
        "eigenvalues_5smallest": evals[:5].tolist(),
        "n_found_beyond": len(found_beyond),
        "n_found_within": len(found_within),
        "matched_tight": matched_tight,
        "matched_loose": matched_loose,
        "n_shown": n_shown,
        "comparisons": comparisons[:60],
    }
    with open(os.path.join(RESULTS_DIR, "connes_exact_weil.json"), "w") as f:
        json.dump(results, f, indent=2)
    print("Results saved.")


if __name__ == "__main__":
    main()
