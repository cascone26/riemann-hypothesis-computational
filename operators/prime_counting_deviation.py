#!/usr/bin/env python3
"""
prime_counting_deviation.py — Connect Connes-τ resonances to prime counting peaks.

Hypothesis: deep resonances in λ_min(τ) correspond to local PEAKS in |ψ(x) - x|.

The mode-alignment condition γ_k × ln(x) ≈ 2πm_k means e^{iγ_k ln(x)} ≈ 1,
which means the zero's contribution to the von Mangoldt explicit formula is
REAL AND POSITIVE: x^{ρ_k}/ρ_k ≈ x^{1/2}/ρ_k (maximally constructive).

Verify by:
1. Computing ψ(x) = Σ_{pk≤x} Λ(k) exactly using sympy/mpmath prime counting
2. Comparing |ψ(x) - x| at resonant vs non-resonant x in range 560-1000

If |ψ(x)-x| is locally elevated at resonant x values, the connection is confirmed.
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 30
OUTFILE = '/Users/scones/projects/riemann/results/prime_counting_deviation.out'

# Deep resonances from scan (lam0 near 0, dt>30s)
RESONANT = {625, 635, 665, 695, 715, 780, 805, 870, 895, 900, 960}

def von_mangoldt_psi(x):
    """
    ψ(x) = Σ_{k≤x, k is prime power} log(prime_base(k))
    Exact computation using integer arithmetic.
    """
    x = int(x)
    psi = 0.0
    for k in range(2, x + 1):
        n, p, f, ok = k, 0, 2, True
        while f * f <= n:
            if n % f == 0:
                if p == 0: p = f
                elif p != f: ok = False; break
                while n % f == 0: n //= f
            f += 1
        if n > 1:
            if p == 0: p = n
            elif p != n: ok = False
        if ok and p > 0:
            psi += float(mp.log(p))
    return psi

def scan():
    """Compute |ψ(x) - x| at resonant and nearby non-resonant x."""
    # Use only step-5 points matching our lam_min scan
    x_vals = list(range(560, 1001, 5))

    print("Connecting deep resonances to prime counting deviation")
    print(f"{'x':>5} {'psi(x)':>12} {'x':>10} {'dev=psi-x':>14} {'|dev|/sqrt(x)':>14} {'deep?'}")
    print("-" * 75)

    results = []
    for x in x_vals:
        psi_x = von_mangoldt_psi(x)
        dev = psi_x - x
        normalized_dev = abs(dev) / np.sqrt(x)
        is_resonant = x in RESONANT
        flag = " ← DEEP" if is_resonant else ""
        line = f"{x:5d} {psi_x:12.2f} {float(x):10.2f} {dev:14.6f} {normalized_dev:14.6f}{flag}"
        if is_resonant or abs(dev) > 20:
            print(line)
        results.append((x, psi_x, dev, normalized_dev, is_resonant))

    # Statistics
    dev_resonant = [abs(r[2]) for r in results if r[4]]
    dev_nonresonant = [abs(r[2]) for r in results if not r[4]]
    norm_resonant = [r[3] for r in results if r[4]]
    norm_nonresonant = [r[3] for r in results if not r[4]]

    print()
    print(f"Mean |ψ(x)-x| at resonant x: {np.mean(dev_resonant):.4f}")
    print(f"Mean |ψ(x)-x| at non-resonant x: {np.mean(dev_nonresonant):.4f}")
    print(f"Mean |ψ(x)-x|/√x at resonant x: {np.mean(norm_resonant):.4f}")
    print(f"Mean |ψ(x)-x|/√x at non-resonant x: {np.mean(norm_nonresonant):.4f}")
    print()
    print("TOP 20 by |ψ(x)-x| (any x in range):")
    results_sorted = sorted(results, key=lambda r: -abs(r[2]))
    for x, psi, dev, norm, is_res in results_sorted[:20]:
        flag = "DEEP" if is_res else "----"
        print(f"  x={x:5d}: |dev|={abs(dev):8.3f}, norm={norm:.4f} [{flag}]")

    with open(OUTFILE, 'w') as f:
        f.write("prime_counting_deviation.py results\n\n")
        f.write(f"Mean |ψ(x)-x| resonant: {np.mean(dev_resonant):.4f}\n")
        f.write(f"Mean |ψ(x)-x| non-resonant: {np.mean(dev_nonresonant):.4f}\n")
        f.write(f"Mean |ψ(x)-x|/√x resonant: {np.mean(norm_resonant):.4f}\n")
        f.write(f"Mean |ψ(x)-x|/√x non-resonant: {np.mean(norm_nonresonant):.4f}\n\n")
        for x, psi, dev, norm, is_res in results:
            flag = "DEEP" if is_res else "----"
            f.write(f"x={x:5d}: psi={psi:.2f} dev={dev:.4f} norm={norm:.4f} [{flag}]\n")
    print(f"\nResults written to {OUTFILE}")

if __name__ == "__main__":
    scan()
