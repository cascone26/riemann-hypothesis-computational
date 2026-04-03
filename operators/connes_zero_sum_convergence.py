"""
Weil Zero-Sum Convergence Test

By the Weil Explicit Formula, the exact Weil matrix τ is:
    τ(V_n, V_m) = Σ_{ALL γ} ê_n(γ) · ê_m(γ)*

where the sum is over ALL non-trivial zeros of ζ.

This means: no separate boundary/prime terms needed — the zero-sum
IS the full Weil form (they cancel perfectly when all zeros are included).

As we add more zeros (K = 10, 50, 100, 500, 1000, 5000, 10000):
  - tau_K is always positive semidefinite (sum of rank-1 PSD matrices)
  - epsilon_K = min eigenvalue of tau_K ≥ 0 always
  - tau_K → tau (true Weil form) as K → ∞
  - epsilon_K → epsilon_∞ = min eigenvalue of true tau

If epsilon_∞ = 0 ↔ RH true, we test: does epsilon_K → 0 as K → ∞?

Key: this approach requires NO approximation of the archimedean term.
The price: we need MANY zeros to converge.

NOTE: With N=80 (161 basis functions), the zero-sum with K≥161 zeros
can have minimum eigenvalue > 0 if the ê vectors don't span everything.
"""

import numpy as np
import mpmath
import json
import os
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR    = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
PLOTS_DIR   = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")

N_TRUNC   = 80
LAMBDA_SQ = 14
L         = float(2 * mpmath.log(mpmath.sqrt(LAMBDA_SQ)))


def load_zeros(count=15000):
    """Load from 2M file (high count), fallback to 100K."""
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
        if len(zeros) >= min(count, 1000):
            print(f"Loaded {len(zeros)} zeros from 2M file (max T={zeros[-1]:.1f})")
            return np.array(zeros)
    except Exception as e:
        print(f"2M file failed: {e}, using 100K")
    path = os.path.join(DATA_DIR, "zeros_100k.txt")
    zeros = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i >= count:
                break
            line = line.strip()
            if line:
                zeros.append(float(line))
    print(f"Loaded {len(zeros)} zeros from 100K file (max T={zeros[-1]:.1f})")
    return np.array(zeros)


def build_zero_sum(gammas, L=L, N=N_TRUNC):
    """
    Build tau = Σ_γ ê_n(γ) ê_m*(γ) vectorized.
    Returns (dim, dim) real symmetric matrix (positive semidefinite).
    """
    ns           = np.arange(-N, N + 1)        # shape (dim,)
    freq_matrix  = gammas[:, None] + 2*np.pi*ns[None, :] / L   # (K, dim)
    half_arg     = freq_matrix * L / 2.0
    safe_freq    = np.where(np.abs(freq_matrix) < 1e-12, 1.0, freq_matrix)
    e_hat        = ((1/np.sqrt(L))
                    * np.exp(1j * half_arg)
                    * 2 * np.sin(half_arg)
                    / safe_freq)
    e_hat[np.abs(freq_matrix) < 1e-12] = np.sqrt(L)  # L'Hopital limit

    # tau = e_hat^T @ conj(e_hat) — shape (dim, dim)
    tau = np.real(e_hat.T @ e_hat.conj())
    return tau


def main():
    print("WEIL ZERO-SUM CONVERGENCE TEST")
    print("=" * 60)
    print(f"N={N_TRUNC}, λ²={LAMBDA_SQ}, L={L:.4f}")
    print("τ = Σ_γ ê_n(γ)ê_m*(γ)  [NO boundary/prime terms]")
    print()

    all_zeros = load_zeros(15000)
    print()

    K_values = [50, 100, 200, 500, 1000, 2000, 5000, 10000]
    K_values = [k for k in K_values if k <= len(all_zeros)]

    results = []
    print(f"{'K':>7}  {'T_max':>8}  {'ε_K':>12}  {'5 smallest eigenvalues':>40}  {'time':>5}")
    print("-" * 85)

    for K in K_values:
        t0      = time.time()
        gammas  = all_zeros[:K]
        tau     = build_zero_sum(gammas)
        evals   = np.linalg.eigvalsh(tau)
        eps     = float(evals[0])
        t       = time.time() - t0

        evals5_str = "  ".join(f"{e:.4f}" for e in evals[:5])
        print(f"{K:>7}  {gammas[-1]:>8.1f}  {eps:>+12.6f}  {evals5_str:>40}  {t:.1f}s")

        results.append({
            "K": K, "T_max": float(gammas[-1]),
            "epsilon_K": eps,
            "5_smallest": evals[:5].tolist(),
        })

    print()
    print("=" * 60)
    print("ANALYSIS")
    print("=" * 60)

    epss = [r["epsilon_K"] for r in results]
    Ks   = [r["K"]         for r in results]

    # Is ε_K decreasing?
    if epss[-1] < epss[0]:
        trend = "DECREASING → 0 (consistent with RH)"
    elif abs(epss[-1] - epss[0]) < 0.01:
        trend = "FLAT (not converging to 0 with these zeros)"
    else:
        trend = "INCREASING (unexpected)"

    print(f"Trend: {trend}")
    print(f"ε at K={Ks[0]}: {epss[0]:.6f}")
    print(f"ε at K={Ks[-1]}: {epss[-1]:.6f}")
    print(f"Change: {epss[-1]-epss[0]:+.6f}")

    # Fit log decay: ε_K ~ A / log(K)^β or ε_K ~ A * K^β
    log_Ks   = np.log(Ks)
    log_epss = np.log(np.abs(epss) + 1e-10)
    try:
        coeffs = np.polyfit(log_Ks, log_epss, 1)
        beta   = coeffs[0]
        A      = np.exp(coeffs[1])
        print(f"Power-law fit: ε_K ≈ {A:.4f} · K^{beta:.4f}")
        # Extrapolate to K=100000, 1000000
        for K_ext in [100000, 1000000, 1e9]:
            print(f"  ε extrapolated at K={K_ext:.0e}: {A * K_ext**beta:.6f}")
    except Exception:
        pass

    print()
    # Key comparison: ε_K from zero-sum vs ε_N from direct Weil (our ~0.32)
    print(f"Direct Weil matrix (approximated archimedean): ε_N ≈ 0.325 (N=150)")
    print(f"Zero-sum at K={Ks[-1]}: ε_K = {epss[-1]:.6f}")
    print()
    if epss[-1] < 0.1:
        print("→ Zero-sum converges fast — direct Weil approx is too pessimistic")
    elif epss[-1] < 0.32:
        print("→ Zero-sum lower than direct Weil — zero-sum is more accurate")
    else:
        print("→ Similar values — both approaches converging to same limit")

    # ── Plot ────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.semilogx(Ks, epss, 'b^-', markersize=8, linewidth=2, label='ε_K (zero-sum)')
    ax.axhline(0.325, color='orange', linestyle='--', linewidth=1.5,
               label='Direct Weil ε_N ≈ 0.325 (N=150)')
    ax.axhline(0, color='red', linestyle='--', linewidth=1.5, label='ε = 0 (RH target)')
    ax.set_xlabel('K (number of zeros in sum)')
    ax.set_ylabel('ε_K (min eigenvalue of zero-sum)')
    ax.set_title(f'Zero-Sum τ Convergence (N={N_TRUNC}, λ²={LAMBDA_SQ})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.loglog(Ks, np.abs(epss) + 1e-10, 'b^-', markersize=8, linewidth=2)
    ax.axhline(0.01, color='green', linestyle=':', linewidth=1, label='0.01 threshold')
    ax.set_xlabel('K (log scale)')
    ax.set_ylabel('|ε_K| (log scale)')
    ax.set_title('Log-Log: Power Law in K?')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.suptitle(
        f'Weil Zero-Sum: τ = Σ_γ ê_n ê_m* | {trend}',
        fontsize=10
    )
    plt.tight_layout()
    out = os.path.join(PLOTS_DIR, "connes_zero_sum_convergence.png")
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"Plot saved: {out}")

    with open(os.path.join(RESULTS_DIR, "connes_zero_sum_convergence.json"), "w") as f:
        json.dump({"N": N_TRUNC, "lambda_sq": LAMBDA_SQ, "L": L,
                   "results": results, "trend": trend}, f, indent=2)
    print("Results saved.")


if __name__ == "__main__":
    main()
