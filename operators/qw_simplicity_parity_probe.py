#!/usr/bin/env python3
"""
qw_simplicity_parity_probe.py  --  Connes 2026 §1174 LEAD probe.

Goal: Check numerically whether the smallest eigenvalue ε(λ) of the Weil
quadratic form QW_λ is SIMPLE (nondegenerate) with an EVEN eigenvector, at
increasing dimension and λ.

This directly addresses Connes' remaining step (§1174): show that the
smallest eigenvalue is simple and even.

Key metrics:
  - Simplicity: gap to 2nd eigenvalue ε₂(λ) - ε₁(λ)
  - Parity: ||even_projected - original|| after even projection
  - Reports for increasing N and λ values

REPORT tier (numerical support); proof requires prolate transfer per Connes.
"""

import sys, json, time
import mpmath as mp
import numpy as np

mp.mp.dps = 40

# Test across multiple lambda and N values
LAMBDA_SQS = [4, 6, 8, 10, 14, 20, 30]  # λ² values (λ in [2, √30] approx)
N_RANGE = [10, 15, 20, 25, 30, 35, 40]   # dimension parameter

# ---------------------------------------------------------------------------
# von Mangoldt data for k <= e^L = lambda^2
# ---------------------------------------------------------------------------
def vonmangoldt_table(kmax):
    tab = []
    for k in range(2, int(kmax) + 1):
        n = k
        p = 0
        f = 2
        ok = True
        while f * f <= n:
            if n % f == 0:
                if p == 0:
                    p = f
                elif p != f:
                    ok = False
                    break
                while n % f == 0:
                    n //= f
            f += 1
        if n > 1:
            if p == 0:
                p = n
            elif p != n:
                ok = False
        if ok and p > 0:
            tab.append((k, mp.log(p)))   # (k, Lambda(k))
    return tab

# ---------------------------------------------------------------------------
# Archimedean density (zero-free)
# ---------------------------------------------------------------------------
def rho(y):
    return mp.e ** (y / 2) / (mp.e ** y - mp.e ** (-y))

def rho_reg(y):
    if y == 0:
        return mp.mpf(1) / 4
    return rho(y) - 1 / (2 * y)

# ---------------------------------------------------------------------------
# Transforms
# ---------------------------------------------------------------------------
def build_tau_fast(lambda_sq, N):
    """Build Weil matrix tau for given λ² and N."""
    L = mp.log(lambda_sq)
    TWO_PI = 2 * mp.pi
    VM = vonmangoldt_table(lambda_sq)

    def _w(k):
        return TWO_PI * k / L

    def Ssin(k):
        if k == 0:
            return mp.mpf(0)
        w = _w(k)
        f = lambda y: mp.sin(w * y) * rho(y)
        return mp.quad(f, [0, L / 2, L])

    def Scos_diag(k):
        w = _w(k)
        g = lambda y: (1 - y / L) * mp.cos(w * y)
        reg_part = mp.quad(lambda y: g(y) * rho_reg(y), [0, L / 2, L])
        head_smooth2 = mp.quad(lambda y: (g(y) - 1) / (2 * y) if y != 0 else
                               ( -1/L ), [0, L/2, L])
        C_R = (mp.log(4 * mp.pi) + mp.euler) / 2
        return reg_part + head_smooth2 + C_R

    def Psin(k):
        s = mp.mpf(0)
        w = _w(k)
        for kk, Lam in VM:
            y = mp.log(kk)
            s += Lam * kk ** (mp.mpf(-1) / 2) * mp.sin(w * y)
        return s

    def Pcos_diag(k):
        s = mp.mpf(0)
        w = _w(k)
        for kk, Lam in VM:
            y = mp.log(kk)
            s += Lam * kk ** (mp.mpf(-1) / 2) * 2 * (1 - y / L) * mp.cos(w * y)
        return s

    def W02(n, m):
        num = (L ** 2 - 16 * mp.pi ** 2 * m * n)
        den = (L ** 2 + 16 * mp.pi ** 2 * m ** 2) * (L ** 2 + 16 * mp.pi ** 2 * n ** 2)
        return 32 * L * mp.sinh(L / 4) ** 2 * num / den

    idx = list(range(-N, N + 1))
    D = 2 * N + 1

    # generators
    b = {}
    a = {}
    for k in range(0, 2 * N + 1):
        bk = -(Ssin(k) - Psin(k)) / mp.pi
        b[k] = bk
        b[-k] = -bk
    for k in range(0, N + 1):
        ak = (Scos_diag(k) - Pcos_diag(k)) + W02(k, k)
        a[k] = ak
        a[-k] = ak

    tau = mp.zeros(D, D)
    for ii, n in enumerate(idx):
        for jj, m in enumerate(idx):
            if n == m:
                tau[ii, jj] = a[n]
            else:
                tau[ii, jj] = (b[n] - b[m]) / (n - m) + W02(n, m)

    # symmetrize
    for ii in range(D):
        for jj in range(ii + 1, D):
            v = (tau[ii, jj] + tau[jj, ii]) / 2
            tau[ii, jj] = v
            tau[jj, ii] = v

    return tau, idx

# ---------------------------------------------------------------------------
# Simplicity and parity checks
# ---------------------------------------------------------------------------
def analyze_smallest_eigenvalue(lambda_sq, N):
    """
    Build QW_λ (tau), compute eigendecomposition, and analyze:
      - Simplicity: gap between 1st and 2nd smallest eigenvalue
      - Parity: how even is the eigenvector?

    Returns: (ε₁, ε₂, gap, eigenvector, parity_score)
    """
    t0 = time.time()
    tau, idx = build_tau_fast(lambda_sq, N)
    D = 2 * N + 1

    # Convert to float array and ensure symmetry
    A = np.array([[float(tau[i, j]) for j in range(D)] for i in range(D)],
                 dtype=float)
    A = (A + A.T) / 2

    # Eigendecomposition (all eigenvalues, sorted)
    evals, evecs = np.linalg.eigh(A)

    # Two smallest eigenvalues
    ε1 = float(evals[0])
    ε2 = float(evals[1]) if len(evals) > 1 else ε1 + 1.0
    gap = ε2 - ε1

    # Eigenvector of smallest eigenvalue
    v = evecs[:, 0]

    # Parity analysis: project to even part
    # even_j = (v_j + v_{-j})/2
    even = np.array([(v[i] + v[D - 1 - i]) / 2 for i in range(D)])

    # Parity score: how much was lost in projection?
    # If v is even, even ≈ v_projected_to_even
    # Compute ||v - even_projected||/||v||
    parity_loss = np.linalg.norm(v - even) / np.linalg.norm(v) if np.linalg.norm(v) > 0 else 1.0

    # More refined: check if v is already even
    # If v_j ≈ v_{D-1-j}, then parity loss should be ~0
    is_even_check = max([abs(v[i] - v[D - 1 - i]) for i in range(D)])

    t_elapsed = time.time() - t0

    return {
        "lambda_sq": lambda_sq,
        "N": N,
        "D": D,
        "epsilon_1": ε1,
        "epsilon_2": ε2,
        "gap": gap,
        "gap_relative": gap / abs(ε1) if ε1 != 0 else 1.0,
        "parity_loss_norm": parity_loss,
        "parity_is_even_max_diff": is_even_check,
        "time_seconds": t_elapsed
    }

# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------
def run_sweep():
    results = {
        "metadata": {
            "title": "Connes 2026 §1174: Smallest eigenvalue simplicity & parity",
            "target": "Is ε(λ) simple (gap > 0) with even eigenvector?",
            "tier": "REPORT (numerical support; proof needs prolate transfer)",
            "timestamp": time.time()
        },
        "runs": []
    }

    print("\n" + "="*80)
    print("QW_λ SIMPLICITY & PARITY PROBE (Connes §1174)")
    print("="*80)
    print(f"{'λ²':>6} {'N':>3} {'D':>4} {'ε₁':>16} {'ε₂':>16} {'gap':>16} {'gap/ε₁':>12} {'parity_loss':>12} {'time':>8}")
    print("-"*80)

    for lambda_sq in LAMBDA_SQS:
        for N in N_RANGE:
            try:
                result = analyze_smallest_eigenvalue(lambda_sq, N)
                results["runs"].append(result)

                # Print row
                lsq = f"{lambda_sq}"
                n_str = f"{N}"
                d_str = f"{result['D']}"
                eps1_str = f"{result['epsilon_1']:.6e}"
                eps2_str = f"{result['epsilon_2']:.6e}"
                gap_str = f"{result['gap']:.6e}"
                gap_rel = f"{result['gap_relative']:.4e}"
                parity = f"{result['parity_loss_norm']:.4e}"
                time_str = f"{result['time_seconds']:.1f}s"

                print(f"{lsq:>6} {n_str:>3} {d_str:>4} {eps1_str:>16} {eps2_str:>16} "
                      f"{gap_str:>16} {gap_rel:>12} {parity:>12} {time_str:>8}")
            except Exception as e:
                print(f"{lambda_sq:>6} {N:>3}  ERROR: {str(e)[:50]}")
                continue

    print("\n" + "="*80)

    # Summary analysis
    print("\nSUMMARY STATISTICS:")
    print("-"*80)

    # Check simplicity: are gaps consistently > 0 and bounded away from 0?
    gaps = [r["gap"] for r in results["runs"]]
    print(f"Simplicity (gap ε₂ - ε₁):")
    print(f"  Min gap: {min(gaps):.6e}")
    print(f"  Max gap: {max(gaps):.6e}")
    print(f"  Mean gap: {np.mean(gaps):.6e}")
    print(f"  All gaps > 0? {all(g > 0 for g in gaps)}")

    # Check parity: how well do the eigenvectors project to even?
    parities = [r["parity_loss_norm"] for r in results["runs"]]
    print(f"\nParity (eigenvector evenness):")
    print(f"  Min parity loss: {min(parities):.6e}")
    print(f"  Max parity loss: {max(parities):.6e}")
    print(f"  Mean parity loss: {np.mean(parities):.6e}")
    print(f"  Median parity loss: {np.median(parities):.6e}")

    # Scaling: does gap grow or shrink with N?
    by_lambda = {}
    for r in results["runs"]:
        lsq = r["lambda_sq"]
        if lsq not in by_lambda:
            by_lambda[lsq] = []
        by_lambda[lsq].append(r)

    print(f"\nScaling with N (for each λ²):")
    for lsq in sorted(by_lambda.keys()):
        runs = sorted(by_lambda[lsq], key=lambda r: r["N"])
        gaps_vs_N = [(r["N"], r["gap"]) for r in runs]
        print(f"  λ² = {lsq:>3}: gaps = {[f'{g:.2e}' for n,g in gaps_vs_N[:3]]} ... "
              f"{[f'{g:.2e}' for n,g in gaps_vs_N[-2:]]}")

    print("\n" + "="*80)
    print("INTERPRETATION (REPORT-TIER):")
    print("-"*80)
    print("• Simplicity: ε₁(λ) is SIMPLE if ε₂ - ε₁ is bounded away from 0")
    print("  → Gap > 1e-10 reliably: SIMPLE ✓")
    print("  → Gap → 0 as N→∞ or λ→∞: NOT simple (degenerate) ✗")
    print()
    print("• Parity: eigenvector is EVEN if parity_loss ≈ 0")
    print("  → Loss < 1e-6: eigenvector IS even ✓")
    print("  → Loss > 0.1: eigenvector NOT even ✗")
    print()
    print("• This is REPORT-TIER numerical evidence.")
    print("• Proof requires HANDLE-TIER: analytical proof via prolate transfer")
    print("  (as Connes §1174 notes: 'analogue is known for prolate operator')")
    print("="*80 + "\n")

    return results

# ---------------------------------------------------------------------------
if __name__ == "__main__":
    results = run_sweep()

    # Save results
    op = "/Users/scones/projects/riemann/results/qw_simplicity_parity.json"
    with open(op, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"Results saved to: {op}\n")
