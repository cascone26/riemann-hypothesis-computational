#!/usr/bin/env python3
"""
zero_accuracy_study.py — Study the digit accuracy pattern across k for the CvS zeros.

From c=13 T=800 and c=23 T=800:
  - k=9 (gamma9=48.005) is consistently bad (35 digits vs 55-80 for others)
  - Odd-k zeros tend to have slightly better accuracy at c=23

This studies:
1. What is sin(gamma_k * L/2) for each zero? (F_cauchy has sin(zL/2) factor)
2. For which k is sin(gamma_k * L/2) ≈ 0? These zeros would be hardest to extract.
3. Does the digit pattern match 1/|sin(gamma_k * L/2)|?

Also: what does the accuracy pattern predict for c=100 T=800 dps=150?
"""
import json
import math
import mpmath as mp

mp.mp.dps = 200

# Known Riemann zeros (high precision)
ZEROS = [
    mp.mpf("14.134725141734693790457251983562470270784257115699"),
    mp.mpf("21.022039638771554992628479593896902777334340524903"),
    mp.mpf("25.010857580145688763213790992562821818659549672558"),
    mp.mpf("30.424876125859513210311897530584091320181560023715"),
    mp.mpf("32.935061587739189690662368964074903488812715603517"),
    mp.mpf("37.586178158825671257217763480705332821405597350831"),
    mp.mpf("40.918719012147495187398126914633254395726165962777"),
    mp.mpf("43.327073280914999519496122165406805782645668371837"),
    mp.mpf("48.005150881167159727942472749427516588896474580065"),
    mp.mpf("49.773832477672302181916784678563724057723178299677"),
]

ZERO_NAMES = [f"γ{k+1}={float(z):.3f}" for k, z in enumerate(ZEROS)]


def analyze_c(c, T_used, results_json=None):
    L = mp.log(c)
    print(f"\n=== c={c}, T={T_used}, L={float(L):.4f} ===")
    print(f"{'k':>3}  {'gamma_k':>10}  {'sin(g*L/2)':>14}  {'|sin|':>10}  {'log10|sin|':>12}", end='')
    if results_json:
        print(f"  {'digits':>8}  {'err':>12}", end='')
    print()

    for k, gamma in enumerate(ZEROS, 1):
        gL2 = gamma * L / 2
        s = mp.sin(gL2)
        log_s = float(mp.log10(abs(s))) if s != 0 else -100
        row = f"{k:>3}  {float(gamma):>10.5f}  {float(s):>14.8f}  {float(abs(s)):>10.6f}  {log_s:>12.4f}"
        if results_json:
            z = results_json['zeros'][k-1] if k <= len(results_json['zeros']) else None
            if z:
                digits = z.get('digits', float('nan'))
                err = z.get('error', 'N/A')
                row += f"  {digits:>8.1f}  {err:>12}"
        print(row)

    # Key: which zeros have |sin(gamma_k * L/2)| very small?
    print(f"\nSmallest |sin(gamma_k * L/2)|:")
    small_sin = sorted(enumerate(ZEROS, 1), key=lambda kz: abs(mp.sin(kz[1] * L / 2)))
    for k, gamma in small_sin[:5]:
        s = mp.sin(gamma * L / 2)
        print(f"  k={k} gamma={float(gamma):.3f}: |sin|={float(abs(s)):.6f}")


def main():
    # Load our results
    c13_data = c23_data = None
    try:
        with open('../results/cvs_highprec_c13.json') as f:
            c13_data = json.load(f)
    except:
        pass
    try:
        with open('../results/cvs_highprec_c23.json') as f:
            c23_data = json.load(f)
    except:
        pass

    analyze_c(13, 800, c13_data)
    analyze_c(23, 800, c23_data)

    # Predict for c=100
    print("\n=== PREDICTION: c=100, T=800 ===")
    c = 100
    L = mp.log(c)
    print(f"L = log(100) = {float(L):.4f}")
    print(f"\nExpected |sin(gamma_k * L/2)| at c=100:")
    svals = []
    for k, gamma in enumerate(ZEROS, 1):
        s = mp.sin(gamma * L / 2)
        svals.append((k, float(gamma), float(abs(s))))
        print(f"  k={k} gamma={float(gamma):.3f}: |sin|={float(abs(s)):.6f}")

    # Zeros with smallest |sin| at c=100 will be hardest
    sorted_by_sin = sorted(svals, key=lambda x: x[2])
    print(f"\nHardest zeros to extract at c=100 (smallest |sin(gamma*L/2)|):")
    for k, gamma, sv in sorted_by_sin[:3]:
        print(f"  k={k} gamma={gamma:.3f}: |sin|={sv:.6f}")


if __name__ == "__main__":
    import os
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    main()
