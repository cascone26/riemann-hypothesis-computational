#!/usr/bin/env python3
"""
multifractal_spectrum.py — PRISM EXTREME item #2: Full multifractal f(α) spectrum.

Extends the single heat-kernel dimension d≈2.46 to the full multifractal f(α) spectrum.
Uses existing heat_kernel_analysis.json data (computed at LAMBDA_SQ=100k and 2M on Olympus).

Two approaches:
A. SCALE-DEPENDENT DIMENSION (from existing JSON data):
   The d(t) variation across scales IS the multifractal signal.
   Maps to the f(α) curve via the Legendre formalism.

B. RÉNYI DIMENSION D_q SWEEP:
   Compute D_q = lim_{ε→0} (1/(q-1)) * log(Σ_i p_i^q) / log(ε)
   where p_i = fraction of spectral measure in box i.
   τ(q) = (q-1)*D_q, f(α) = Legendre transform.
   Uses dim_evolution_scan.json: lam_min(lambda_sq) data as the spectral measure.

A monofractal has d(t) = const and f(α) is a point.
A multifractal has d(t) varying and f(α) is a curved parabola-like shape.
"""

import numpy as np
import json
import time

def load_heat_kernel():
    """Load scale-dependent dimension data from existing computation."""
    d = json.load(open("/Users/scones/projects/riemann/results/heat_kernel_analysis.json"))
    results = {}
    for key in ['100k', '2M']:
        if key in d:
            entry = d[key]
            scale_d = entry.get('scale_dependent_d', [])
            results[key] = {
                'd_simple': entry['d_simple'],
                'scale_d': scale_d  # [[t_lo, t_hi, d], ...]
            }
    return results

def legendre_transform_from_dq(q_vals, tau_q):
    """
    Given τ(q) = (q-1)*D_q, compute f(α) spectrum.
    α(q) = -dτ/dq, f(α) = q*α - τ(q).
    """
    alpha_q = -np.gradient(tau_q, q_vals)  # NOTE: sign depends on convention
    # Standard convention: τ(q) from Z_q ~ ε^τ(q), α = dτ/dq, f = qα - τ
    alpha_q_std = np.gradient(tau_q, q_vals)
    f_alpha = q_vals * alpha_q_std - tau_q
    return alpha_q_std, f_alpha

def spectral_multifractal_from_scan():
    """
    Use dim_evolution_scan.json to compute Rényi dimensions D_q
    of the spectral measure λ_min(lambda_sq).

    Treat |lam_min(lambda_sq)| as a measure density on the lambda_sq axis.
    Box-count at different scales to get D_q.
    """
    scan = json.load(open("/Users/scones/projects/riemann/results/dim_evolution_scan.json"))

    # Extract lambda_sq → lam_min pairs
    data = [(s['lambda_sq'], abs(s['lam_min'])) for s in scan if s['lam_min'] != 0]
    lsq_vals = np.array([x[0] for x in data])
    lam_vals = np.array([x[1] for x in data])

    # Normalize to probability distribution (spectral measure)
    total = lam_vals.sum()
    p = lam_vals / total

    # Compute generalized dimensions D_q at multiple scales (box sizes)
    q_vals = np.linspace(-4, 4, 33)

    # Use box-counting at 3 different resolutions on [lambda_sq_min, lambda_sq_max]
    lsq_min, lsq_max = lsq_vals.min(), lsq_vals.max()

    n_boxes_list = [4, 8, 16]
    D_q_by_scale = {}

    for n_boxes in n_boxes_list:
        box_edges = np.linspace(lsq_min, lsq_max, n_boxes + 1)
        box_p = np.zeros(n_boxes)
        for b in range(n_boxes):
            mask = (lsq_vals >= box_edges[b]) & (lsq_vals < box_edges[b+1])
            box_p[b] = p[mask].sum()
        box_p = box_p[box_p > 0]  # remove empty boxes
        eps = (lsq_max - lsq_min) / n_boxes

        D_q_row = []
        for q in q_vals:
            if abs(q - 1) < 1e-6:
                # D_1 = Shannon entropy / log(1/eps)
                S = -np.sum(box_p * np.log(box_p + 1e-300))
                D_q_row.append(S / np.log(1/eps))
            else:
                Z_q = np.sum(box_p**q)
                D_q_row.append(np.log(Z_q) / ((q-1) * np.log(eps)))
        D_q_by_scale[n_boxes] = np.array(D_q_row)

    # Average D_q across scales (should be scale-independent if truly multifractal)
    D_q_mean = np.mean([D_q_by_scale[n] for n in n_boxes_list], axis=0)
    tau_q = (q_vals - 1) * D_q_mean

    return q_vals, D_q_mean, tau_q, lsq_vals, lam_vals

if __name__ == "__main__":
    print("Multifractal f(α) spectrum — CvS operator")
    print("==========================================")
    print()

    # === Part A: Scale-dependent d from existing heat kernel data ===
    print("=== Part A: Scale-dependent dimension (existing heat_kernel_analysis.json) ===")
    hk = load_heat_kernel()

    for key, info in hk.items():
        print(f"\nLAMBDA_SQ = {key}:")
        print(f"  d_simple = {info['d_simple']:.4f}")
        print(f"  Scale-dependent d(t) — the multifractal signal:")

        alpha_vals = []
        for t_lo, t_hi, d_loc in info['scale_d']:
            t_mid = np.sqrt(t_lo * t_hi)  # geometric midpoint
            alpha_loc = d_loc / 2.0  # local Hölder exponent α = d/2
            alpha_vals.append((t_mid, d_loc, alpha_loc))
            print(f"    t ∈ [{t_lo:.4f}, {t_hi:.4f}]:  d = {d_loc:.4f}  →  α = {alpha_loc:.4f}")

        # Width of multifractal spectrum
        alphas = [a[2] for a in alpha_vals]
        print(f"  α range: [{min(alphas):.4f}, {max(alphas):.4f}]")
        print(f"  Δα = {max(alphas)-min(alphas):.4f}  "
              f"({'monofractal (Δα≈0)' if max(alphas)-min(alphas) < 0.1 else 'MULTIFRACTAL (Δα>0)'})")

        # Is d(t) monotone or has extremum?
        d_vals = [a[1] for a in alpha_vals]
        if d_vals[-1] > d_vals[0]:
            print(f"  d(t): INCREASING with t (UV smaller than IR) → scale-dependent spectral geometry")
        elif d_vals[-1] < d_vals[0]:
            print(f"  d(t): DECREASING with t → long-range spectral compression")
        else:
            print(f"  d(t): approximately constant → monofractal")

    print()
    print("=== Part B: Rényi dimensions D_q from spectral measure scan ===")
    print("(Using dim_evolution_scan.json: |λ_min|(lambda_sq) as spectral measure)")

    q_vals, D_q, tau_q, lsq_vals, lam_vals = spectral_multifractal_from_scan()
    alpha_q, f_alpha = legendre_transform_from_dq(q_vals, tau_q)

    print(f"\nSpectral measure: {len(lsq_vals)} points, lambda_sq ∈ [{lsq_vals.min():.1f}, {lsq_vals.max():.1f}]")
    print(f"\nRényi dimension profile D_q:")
    print(f"  D_0 = {D_q[q_vals==0.0][0]:.4f}  (Hausdorff dimension of support)")
    D_1_idx = np.argmin(np.abs(q_vals - 1.0))
    D_2_idx = np.argmin(np.abs(q_vals - 2.0))
    print(f"  D_1 = {D_q[D_1_idx]:.4f}  (information dimension)")
    print(f"  D_2 = {D_q[D_2_idx]:.4f}  (correlation dimension)")
    print(f"  D_q range: [{D_q.min():.4f}, {D_q.max():.4f}]")
    print(f"  ΔD = D_0 - D_∞ = {D_q[0] - D_q[-1]:.4f}  (multifractal width)")

    print(f"\nτ(q) and f(α) spectrum:")
    print(f"  {'q':>6}   {'D_q':>8}   {'τ(q)':>8}   {'α(q)':>8}   {'f(α)':>8}")
    print(f"  {'─'*50}")
    for i in range(0, len(q_vals), 4):
        q = q_vals[i]
        print(f"  {q:>6.2f}   {D_q[i]:>8.4f}   {tau_q[i]:>8.4f}   {alpha_q[i]:>8.4f}   {f_alpha[i]:>8.4f}")

    # Physical region
    phys = f_alpha >= -0.2
    if phys.sum() > 1:
        alpha_phys = alpha_q[phys]
        f_phys = f_alpha[phys]
        print(f"\nf(α) spectrum summary (physical region f(α)≥-0.2):")
        print(f"  α range: [{alpha_phys.min():.4f}, {alpha_phys.max():.4f}]")
        print(f"  Δα = {alpha_phys.max()-alpha_phys.min():.4f}")
        print(f"  α at f(α)_max: {alpha_q[np.argmax(f_alpha)]:.4f}")
        print(f"  f(α)_max = {f_alpha.max():.4f}")

    # Is it multifractal? D_0 > D_1 > D_2 is the signature
    D_0 = D_q[np.argmin(np.abs(q_vals))]
    print(f"\nMultifractality test (should have D_0 > D_1 > D_2 for genuine multifractal):")
    print(f"  D_0 = {D_0:.4f}  {'>' if D_0 > D_q[D_1_idx] else '≤'}  D_1 = {D_q[D_1_idx]:.4f}  "
          f"{'>' if D_q[D_1_idx] > D_q[D_2_idx] else '≤'}  D_2 = {D_q[D_2_idx]:.4f}")
    if D_0 > D_q[D_1_idx] > D_q[D_2_idx]:
        print(f"  → MULTIFRACTAL spectrum confirmed (D_0 > D_1 > D_2)")
    else:
        print(f"  → Simple or monofractal spectrum")

    # Save
    result = {
        "heat_kernel_scale_d": {k: v['scale_d'] for k, v in hk.items()},
        "q_vals": q_vals.tolist(),
        "D_q": D_q.tolist(),
        "tau_q": tau_q.tolist(),
        "alpha_q": alpha_q.tolist(),
        "f_alpha": f_alpha.tolist(),
        "summary": {
            "D_0": float(D_0),
            "D_1": float(D_q[D_1_idx]),
            "D_2": float(D_q[D_2_idx]),
            "delta_D": float(D_q[0] - D_q[-1]),
            "alpha_range_partA_100k": {
                "min": float(min(a[2] for a in alpha_vals)) if hk.get('100k') else None,
                "max": float(max(a[2] for a in alpha_vals)) if hk.get('100k') else None,
            }
        }
    }
    out = "/Users/scones/projects/riemann/results/multifractal_spectrum.json"
    with open(out, "w") as fh:
        json.dump(result, fh, indent=2)
    print(f"\nSaved → {out}")
