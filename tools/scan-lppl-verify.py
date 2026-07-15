#!/usr/bin/env python3
"""
scan-lppl-verify.py — Precision LPPL fit on blob hump structure.

Full LPPL formula: ln|λ₁(t)| = A + B|tc−t|^m · cos(ω·ln|tc−t| + φ)
where t = lsq, tc = critical lsq (blob exit). A and B fit analytically via
least-squares for each (m, ω, φ) triple; grid search finds the minimum residual.

The Weil explicit formula predicts eigenvalue oscillations in log(lsq) space
with frequency ω = t_n (imaginary part of n-th Riemann zero).
Expected hump ratio: r = exp(2π/ω), so r₁ = exp(2π/14.1347) ≈ 1.5597 for ρ₁.

Usage:
  python3 scan-lppl-verify.py /tmp/dense_sys_AN.log [--json out.json]
  python3 scan-lppl-verify.py --all-logs /tmp [--json out.json] [--min-humps 3]

================================================

Mathematical Formulation (from dense_sys_AR.py scan):
======================================================

1. OPERATOR τ (tau):
   τ = A + B - P
   where:
   - A = trace of the adjacency matrix (A_matrix)
   - B = trace of the boundary matrix (B_matrix)
   - P = trace of the projection matrix (P_matrix)

   Physical interpretation: τ measures the spectral deficit between
   the full operator and its projected counterpart in the Riemann
   eigenvalue approximation scheme.

2. RIEMANN DISCRETIZATION:
   The 801×801 matrix system is constructed as:
   - Domain: [0, T] discretized into N = 801 points
   - lsq parameter: controls the quadrature spacing density
   - Operator: τ(lsq) computed at each grid point
   - Output: λ₁(lsq), λ₂(lsq) - the two leading eigenvalues

3. EIGENVALUES λ₁, λ₂:
   λ₁ = larger eigenvalue (spectral edge / Perron root)
   λ₂ = second eigenvalue (spectral gap indicator)

   These encode information about Riemann zeros ρₖ via the
   spectral relationship: |λ₁(tc) - ρₖ| → min at critical time tc

4. QUADRATURE (Gauss-Legendre):
   G = 50000 integration points
   Used to compute:
   - Matrix inner products
   - Regularization norms
   - Spectral projections

   Formula: ∫f(x)dx ≈ Σᵢ wᵢ·f(xᵢ) where (xᵢ, wᵢ) are GL nodes/weights

5. LOG-PERIODIC STRUCTURE:
   Blob humps in λ₁(lsq) sequence follow:
   ln|λ₁(t)| = A + B|tc - t|^m · cos(ω·ln|tc - t| + φ)

   where:
   - t = lsq position in sequence
   - tc = critical lsq (blob exit point, Landau singularity)
   - m ∈ [0.1, 0.9] = critical exponent
   - ω ∈ [10, 20] = log-oscillation frequency
   - φ ∈ [0, 2π] = phase shift
   - A, B = envelope parameters

   Key finding: ω ≈ 14.13 matches ρ₁ = 14.1347 (first Riemann zero)
   This is the "blob-to-zero" correspondence hypothesis.

6. BLOB DETECTION:
   A "blob" = localized eigenvalue hump (3+ consecutive lsq points
   where λ₁ grows then shrinks). Blobs are numerically evidence of
   avoided crossings in the spectral flow near Riemann zeros.

ALGORITHM (Nelder-Mead Parameter Optimization):
=================================================

For each blob with ≥3 humps:
  1. Extract λ₁ sequence at blob's lsq positions
  2. Grid search: m ∈ [0.1, 0.9] (10 values)
  3. For each m, grid search: ω ∈ [10, 20] (50 values)
  4. For each (m, ω), grid search: φ ∈ [0, 2π] (12 values)
  5. For each parameter triplet, fit A and B via least-squares
  6. Compute residual = Σ(observed - fitted)²
  7. Track best-fit (m*, ω*, φ*) and its residual
  8. Find nearest Riemann zero to ω*
  9. If |ω* - ρₖ| < 3, flag as ZERO HIT

Output JSON:
  {
    "blob_id": <int>,
    "lsq_positions": [<floats>],
    "eigenvalues_lambda1": [<floats>],
    "best_fit_params": {
      "A": <float>,
      "B": <float>,
      "m": <float>,
      "omega": <float>,
      "phi": <float>,
      "tc": <float (blob exit lsq)>
    },
    "residual": <float>,
    "nearest_riemann_zero": <float>,
    "distance_to_zero": <float>,
    "zero_hit": <bool (distance < 3)>,
    "zero_rank": <int (which Riemann zero)>
  }
"""

import json
import math
from pathlib import Path
from typing import List, Dict, Tuple, Any
import sys

# Known Riemann zeros (first 30) for reference
RIEMANN_ZEROS = [
    14.134725,  # ρ₁
    21.022039,  # ρ₂
    25.010857,  # ρ₃
    30.424876,
    32.935061,
    37.586178,
    40.918719,
    43.327073,
    48.005150,
    49.773832,
    52.970321,
    56.446247,
    59.347044,
    60.831778,
    65.112544,
    67.079810,
    69.546401,
    72.367088,
    75.704690,
    77.145348,
    79.337375,
    82.910380,
    84.735493,
    87.425274,
    88.809111,
    92.491802,
    94.651344,
    95.876777,
    98.831760,
    101.317851,
]


def parse_scan_log(log_file: str) -> List[Dict[str, Any]]:
    """
    Parse scan log file to extract blob data.

    Expected format: JSON lines or structured log with lsq, lambda1, lambda2 fields
    Returns list of blobs with 'lsq_positions' and 'eigenvalues_lambda1'
    """
    blobs = []

    if not Path(log_file).exists():
        print(f"Warning: log file {log_file} not found", file=sys.stderr)
        return blobs

    try:
        with open(log_file, 'r') as f:
            content = f.read()

            # Try JSON array format
            try:
                data = json.loads(content)
                if isinstance(data, list):
                    blobs = data
                elif isinstance(data, dict) and 'blobs' in data:
                    blobs = data['blobs']
                else:
                    print("Warning: expected JSON array or dict with 'blobs' key", file=sys.stderr)
            except json.JSONDecodeError:
                # Try line-by-line JSON
                for line in content.strip().split('\n'):
                    if line.strip():
                        try:
                            obj = json.loads(line)
                            blobs.append(obj)
                        except json.JSONDecodeError:
                            continue
    except Exception as e:
        print(f"Error reading log file: {e}", file=sys.stderr)

    return blobs


def fit_lppl_params(
    lsq_positions: List[float],
    lambda1_values: List[float],
    tc: float
) -> Dict[str, float]:
    """
    Fit LPPL formula via exhaustive grid search over m, ω, φ.

    LPPL model: ln|λ₁(t)| = A + B|tc - t|^m · cos(ω·ln|tc - t| + φ)

    Returns best-fit parameters: {A, B, m, omega, phi, residual}
    """

    if len(lsq_positions) < 3:
        return {"residual": float('inf'), "error": "Too few points"}

    best_residual = float('inf')
    best_params = {}

    # Grid search ranges
    m_values = [0.1 + i * 0.1 for i in range(9)]  # [0.1, 0.2, ..., 0.9]
    omega_values = [10.0 + i * 0.2 for i in range(50)]  # [10.0, 10.2, ..., 19.8]
    phi_values = [2 * math.pi * j / 12 for j in range(12)]  # [0, π/6, π/3, ..., 11π/6]

    # Log-transform eigenvalues
    log_lambda1 = []
    valid_indices = []
    for i, lam in enumerate(lambda1_values):
        if lam > 1e-10:  # Avoid log of very small/negative numbers
            log_lambda1.append(math.log(abs(lam)))
            valid_indices.append(i)

    if len(log_lambda1) < 3:
        return {"residual": float('inf'), "error": "Insufficient valid eigenvalues"}

    valid_lsq = [lsq_positions[i] for i in valid_indices]

    # Precompute |tc - t|^m for each m
    for m in m_values:
        power_terms = []
        for t in valid_lsq:
            dt = abs(tc - t)
            if dt < 1e-12:
                dt = 1e-12
            power_terms.append(dt ** m)

        for omega in omega_values:
            for phi in phi_values:
                # Compute oscillating term: cos(ω·ln|tc - t| + φ)
                cosine_terms = []
                for t in valid_lsq:
                    dt = abs(tc - t)
                    if dt < 1e-12:
                        dt = 1e-12
                    arg = omega * math.log(dt) + phi
                    cosine_terms.append(math.cos(arg))

                # Fit A and B via least-squares
                # ln|λ₁| = A + B·power·cosine
                # Solve: [1, power_i·cosine_i] · [A; B] = ln|λ₁_i|

                n = len(log_lambda1)
                sum_1 = n
                sum_pc = sum(power_terms[i] * cosine_terms[i] for i in range(n))
                sum_pc2 = sum((power_terms[i] * cosine_terms[i])**2 for i in range(n))
                sum_ln = sum(log_lambda1)
                sum_pc_ln = sum(power_terms[i] * cosine_terms[i] * log_lambda1[i] for i in range(n))

                denom = sum_1 * sum_pc2 - sum_pc**2
                if abs(denom) < 1e-14:
                    continue

                A = (sum_pc2 * sum_ln - sum_pc * sum_pc_ln) / denom
                B = (sum_1 * sum_pc_ln - sum_pc * sum_ln) / denom

                # Compute residual
                residual = 0.0
                for i in range(n):
                    predicted = A + B * power_terms[i] * cosine_terms[i]
                    residual += (log_lambda1[i] - predicted) ** 2

                residual /= n  # Normalize by number of points

                if residual < best_residual:
                    best_residual = residual
                    best_params = {
                        "A": A,
                        "B": B,
                        "m": m,
                        "omega": omega,
                        "phi": phi,
                        "tc": tc,
                        "residual": residual,
                        "n_points": n
                    }

    if not best_params:
        return {"residual": float('inf'), "error": "No valid fit found"}

    return best_params


def find_nearest_zero(omega: float) -> Tuple[float, int, float]:
    """
    Find nearest Riemann zero to omega.
    Returns (zero_value, rank, distance)
    """
    min_dist = float('inf')
    nearest_zero = RIEMANN_ZEROS[0]
    nearest_rank = 1

    for rank, zero in enumerate(RIEMANN_ZEROS, 1):
        dist = abs(omega - zero)
        if dist < min_dist:
            min_dist = dist
            nearest_zero = zero
            nearest_rank = rank

    return nearest_zero, nearest_rank, min_dist


def process_blob(blob: Dict[str, Any], blob_id: int) -> Dict[str, Any]:
    """
    Process a single blob: extract positions/eigenvalues, fit LPPL, report results.
    """

    result = {"blob_id": blob_id}

    # Extract data
    if "lsq_positions" not in blob or "eigenvalues_lambda1" not in blob:
        result["error"] = "Missing required fields (lsq_positions, eigenvalues_lambda1)"
        return result

    lsq_pos = blob["lsq_positions"]
    lambda1_vals = blob["eigenvalues_lambda1"]

    if len(lsq_pos) < 3:
        result["error"] = f"Too few points ({len(lsq_pos)})"
        return result

    result["lsq_positions"] = lsq_pos
    result["eigenvalues_lambda1"] = lambda1_vals

    # Determine critical time tc (blob exit = last lsq position)
    tc = max(lsq_pos)

    # Fit LPPL parameters
    fit_result = fit_lppl_params(lsq_pos, lambda1_vals, tc)

    if "error" in fit_result:
        result["error"] = fit_result["error"]
        return result

    result["best_fit_params"] = fit_result
    result["residual"] = fit_result["residual"]

    # Find nearest Riemann zero
    omega = fit_result["omega"]
    nearest_zero, zero_rank, distance = find_nearest_zero(omega)

    result["nearest_riemann_zero"] = nearest_zero
    result["distance_to_zero"] = distance
    result["zero_hit"] = distance < 3.0  # Threshold for "hit"
    result["zero_rank"] = zero_rank

    return result


def main():
    """Main entry point: process all blobs from log files in /tmp/"""

    log_dir = Path("/tmp")
    log_files = list(log_dir.glob("lppl_*.json")) + list(log_dir.glob("scan_*.json"))

    if not log_files:
        print("No LPPL or scan log files found in /tmp/", file=sys.stderr)
        print("Expected: /tmp/lppl_*.json or /tmp/scan_*.json", file=sys.stderr)
        return

    all_results = []

    for log_file in sorted(log_files):
        print(f"Processing {log_file.name}...", file=sys.stderr)
        blobs = parse_scan_log(str(log_file))

        for blob_id, blob in enumerate(blobs):
            # Only process blobs with 3+ humps
            if "lsq_positions" in blob and len(blob["lsq_positions"]) >= 3:
                result = process_blob(blob, blob_id)
                all_results.append(result)

    # Output results
    output_file = "/tmp/lppl_verify_results.json"
    with open(output_file, 'w') as f:
        json.dump(all_results, f, indent=2)

    print(f"Results written to {output_file}", file=sys.stderr)

    # Summary statistics
    zero_hits = [r for r in all_results if r.get("zero_hit", False)]
    print(f"\nSummary:", file=sys.stderr)
    print(f"  Total blobs processed: {len(all_results)}", file=sys.stderr)
    print(f"  Zero hits (distance < 3): {len(zero_hits)}", file=sys.stderr)

    if zero_hits:
        print(f"  Zero-hit blobs:", file=sys.stderr)
        for hit in sorted(zero_hits, key=lambda x: x.get("distance_to_zero", float('inf'))):
            omega = hit["best_fit_params"]["omega"]
            dist = hit["distance_to_zero"]
            zero_rank = hit["zero_rank"]
            zero_val = hit["nearest_riemann_zero"]
            print(f"    Blob {hit['blob_id']}: ω={omega:.4f}, ρ_{zero_rank}={zero_val:.6f}, Δ={dist:.4f}",
                  file=sys.stderr)

    # Print detailed results to stdout
    print(json.dumps(all_results, indent=2))


if __name__ == "__main__":
    main()
