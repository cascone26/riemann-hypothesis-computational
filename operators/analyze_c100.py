#!/usr/bin/env python3
"""
Analyze c=100 CvS computation results from Olympus.

Reads /home/cobo/riemann/results/cvs_highprec_c100.json via SSH,
extracts Riemann zeros with precision analysis, detects outliers,
and appends findings to SESSION_FINDINGS_2026_07_11.md.

Required:
  - mpmath (for high-precision arithmetic)
  - SSH access to olympus (configured)
"""

import subprocess
import json
import sys
import os
from pathlib import Path
from datetime import datetime
import math

try:
    import mpmath
except ImportError:
    print("ERROR: mpmath not installed. Install with: pip install mpmath")
    sys.exit(1)


# ============================================================================
# CONFIGURATION
# ============================================================================

OLYMPUS_HOST = "olympus"
REMOTE_RESULT_FILE = "/home/cobo/riemann/results/cvs_highprec_c100.json"
REMOTE_CHAIN_LOG = "/home/cobo/riemann/results/chain_log.txt"

LOCAL_PROJECT_ROOT = Path(__file__).parent.parent
LOCAL_RESULTS_DIR = LOCAL_PROJECT_ROOT / "results"
LOCAL_FINDINGS_FILE = LOCAL_RESULTS_DIR / "SESSION_FINDINGS_2026_07_11.md"

# Known Riemann zeros (as strings for mpmath precision)
KNOWN_ZEROS = {
    1: "14.134725141734693790457251983562",
    2: "21.022039638771554992628479593896",
    3: "25.010857580145688763213790992562",
    4: "30.424876125859513210311897530584",
    5: "32.935061587739189690662368964074",
    6: "37.586178158825671257217763480705",
    7: "40.918719012147495187398126914633",
    8: "43.327073280914999519496122165406",
    9: "48.005150881167159727942472749427",
    10: "49.773832477672302181916784678563",
}


# ============================================================================
# SSH DATA FETCH
# ============================================================================

def fetch_remote_file(host, remote_path):
    """
    Fetch file from remote machine via SSH.
    Returns file contents as string, or raises exception on failure.
    """
    try:
        result = subprocess.run(
            ["ssh", host, f"cat {remote_path}"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        if result.returncode != 0:
            raise RuntimeError(f"SSH failed: {result.stderr}")
        return result.stdout
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"SSH timeout connecting to {host}")
    except FileNotFoundError:
        raise RuntimeError("ssh command not found; ensure SSH is configured")


def fetch_json_results(host, remote_path):
    """Fetch and parse JSON results file."""
    content = fetch_remote_file(host, remote_path)
    try:
        data = json.loads(content)
        return data
    except json.JSONDecodeError as e:
        raise RuntimeError(f"Invalid JSON in {remote_path}: {e}")


def fetch_chain_log(host, remote_path):
    """Fetch chain log; returns raw text (may not exist)."""
    try:
        return fetch_remote_file(host, remote_path)
    except Exception as e:
        return f"[Chain log unavailable: {e}]"


# ============================================================================
# JSON PARSING
# ============================================================================

def inspect_json_structure(data):
    """Print top-level keys for debugging."""
    print("JSON structure:")
    if isinstance(data, dict):
        for key in sorted(data.keys()):
            val = data[key]
            if isinstance(val, (list, dict)):
                print(f"  {key}: {type(val).__name__} (len={len(val) if isinstance(val, (list, dict)) else 'N/A'})")
            else:
                print(f"  {key}: {type(val).__name__} = {val}")
    else:
        print(f"  Top level is {type(data).__name__}")


def extract_results(data):
    """
    Extract lam_min, zeros, and other metadata from JSON.
    Returns dict with keys: lam_min, zeros, raw_data.
    """
    # Try to infer structure
    lam_min = None
    zeros_data = None

    # Common keys for lam_min
    for key in ["lam_min", "lambda_min", "ground_state_eigenvalue"]:
        if key in data:
            lam_min = data[key]
            break

    # Common keys for zeros
    for key in ["zeros", "gamma_k", "extracted_zeros", "riemann_zeros"]:
        if key in data:
            zeros_data = data[key]
            break

    # If still not found, try list directly
    if zeros_data is None and isinstance(data, dict):
        # Maybe zeros are stored as a list at top level
        for key, val in data.items():
            if isinstance(val, list) and len(val) > 0 and isinstance(val[0], (int, float)):
                zeros_data = val
                break

    if lam_min is None or zeros_data is None:
        raise ValueError("Could not extract lam_min or zeros from JSON. Structure unknown.")

    return {
        "lam_min": lam_min,
        "zeros": zeros_data,
        "raw_data": data,
    }


# ============================================================================
# PRECISION ANALYSIS
# ============================================================================

def setup_mpmath(dps=200):
    """Configure mpmath for high-precision arithmetic."""
    mpmath.mp.dps = dps


def compute_precision_metrics(extracted_zeros):
    """
    Compare extracted zeros against known reference values.
    Returns list of dicts with keys:
      k, gamma_extracted, gamma_ref, error, digit_precision, sign
    """
    metrics = []
    setup_mpmath(dps=200)

    for k in sorted(extracted_zeros.keys()):
        gamma_ext = mpmath.mpf(str(extracted_zeros[k]))

        if k not in KNOWN_ZEROS:
            continue  # Skip if no reference

        gamma_ref = mpmath.mpf(KNOWN_ZEROS[k])
        error = abs(gamma_ext - gamma_ref)

        # digit_precision = -log10(error)
        if error > 0:
            digit_precision = -mpmath.log10(error)
        else:
            digit_precision = mpmath.mpf("inf")

        # Sign of extracted zero
        sign = "+" if gamma_ext > 0 else "-"

        metrics.append({
            "k": k,
            "gamma_extracted": float(gamma_ext),
            "gamma_ref": float(gamma_ref),
            "error": float(error),
            "digit_precision": float(digit_precision),
            "sign": sign,
        })

    return metrics


def identify_outliers(metrics, sigma_threshold=2.0):
    """
    Identify zeros where digit_precision is >sigma_threshold
    below the mean.
    """
    if not metrics:
        return []

    precisions = [m["digit_precision"] for m in metrics]
    mean_precision = sum(precisions) / len(precisions)
    variance = sum((p - mean_precision) ** 2 for p in precisions) / len(precisions)
    std_dev = math.sqrt(variance)

    outliers = []
    for m in metrics:
        z_score = (mean_precision - m["digit_precision"]) / std_dev if std_dev > 0 else 0
        if z_score > sigma_threshold:
            outliers.append({
                **m,
                "z_score": z_score,
                "threshold": sigma_threshold,
            })

    return outliers


def analyze_k7_outlier(extracted_zeros, lam_min_value):
    """
    Special analysis for k=7 zero.
    Compute |sin(γ₇L/2)| where L = log(100) ≈ 4.605...
    """
    setup_mpmath(dps=200)

    if 7 not in extracted_zeros:
        return None

    gamma7 = mpmath.mpf(str(extracted_zeros[7]))
    L = mpmath.log(100)

    # Compute |sin(γ₇L/2)|
    sin_val = mpmath.sin(gamma7 * L / 2)
    abs_sin_val = abs(sin_val)

    return {
        "gamma7": float(gamma7),
        "L": float(L),
        "sin(gamma7 * L / 2)": float(sin_val),
        "abs_sin(gamma7 * L / 2)": float(abs_sin_val),
        "lam_min": lam_min_value,
    }


# ============================================================================
# OUTPUT FORMATTING
# ============================================================================

def format_findings(results_dict, metrics, outliers, k7_analysis, chain_log_snippet):
    """Format findings as markdown section."""
    section = "## §49. CvS High-Precision Analysis (c=100, T=800, dps=150)\n\n"
    section += f"**Date**: {datetime.now().isoformat()}\n\n"

    lam_min = results_dict["lam_min"]
    sign = "+" if lam_min > 0 else "-"
    section += f"### lam_min\n"
    section += f"- **Value**: {lam_min}\n"
    section += f"- **Sign**: {sign}\n\n"

    section += f"### Extracted Riemann Zeros\n"
    section += f"**Count**: {len(metrics)} zeros\n\n"

    # Table of zeros
    section += "| k | γ_extracted | γ_ref | error | digit_precision | sign |\n"
    section += "|---|---|---|---|---|---|\n"
    for m in sorted(metrics, key=lambda x: x["k"]):
        section += (
            f"| {m['k']} | {m['gamma_extracted']:.15g} | {m['gamma_ref']:.15g} | "
            f"{m['error']:.3e} | {m['digit_precision']:.2f} | {m['sign']} |\n"
        )
    section += "\n"

    # Precision summary
    if metrics:
        precisions = [m["digit_precision"] for m in metrics]
        mean_p = sum(precisions) / len(precisions)
        section += f"### Precision Summary\n"
        section += f"- **Mean digit precision**: {mean_p:.2f} digits\n"
        section += f"- **Min digit precision**: {min(precisions):.2f} digits (zero k={metrics[precisions.index(min(precisions))]['k']})\n"
        section += f"- **Max digit precision**: {max(precisions):.2f} digits (zero k={metrics[precisions.index(max(precisions))]['k']})\n\n"

    # Outliers
    if outliers:
        section += f"### Outliers (>2σ below mean precision)\n"
        section += f"**Count**: {len(outliers)} outliers detected.\n\n"
        section += "| k | z_score | digit_precision | error |\n"
        section += "|---|---|---|---|\n"
        for o in sorted(outliers, key=lambda x: x["k"]):
            section += (
                f"| {o['k']} | {o['z_score']:.2f} | {o['digit_precision']:.2f} | "
                f"{o['error']:.3e} |\n"
            )
        section += "\n"
    else:
        section += "### Outliers\nNo significant outliers detected.\n\n"

    # k=7 analysis
    if k7_analysis:
        section += "### k=7 Special Analysis\n"
        section += f"γ₇ = {k7_analysis['gamma7']:.15g}\n"
        section += f"L = log(100) = {k7_analysis['L']:.15g}\n"
        section += f"|sin(γ₇L/2)| = {k7_analysis['abs_sin(gamma7 * L / 2)']:.3e}\n"
        section += f"lam_min = {k7_analysis['lam_min']}\n\n"

    # Chain log check
    section += "### Chain Log Verification\n"
    if "unavailable" in chain_log_snippet.lower():
        section += f"Chain log: Not accessible or doesn't exist yet.\n\n"
    else:
        section += f"Chain log (first 500 chars):\n```\n{chain_log_snippet[:500]}\n```\n\n"

    section += "---\n\n"
    return section


# ============================================================================
# MAIN
# ============================================================================

def main():
    print(f"[{datetime.now().isoformat()}] Starting c=100 analysis...")

    # Step 1: Fetch results from Olympus
    print(f"[1/5] Fetching results from {OLYMPUS_HOST}:{REMOTE_RESULT_FILE}")
    try:
        results_data = fetch_json_results(OLYMPUS_HOST, REMOTE_RESULT_FILE)
    except Exception as e:
        print(f"ERROR fetching results: {e}")
        sys.exit(1)

    # Step 2: Inspect and extract
    print("[2/5] Inspecting JSON structure...")
    inspect_json_structure(results_data)

    try:
        results_dict = extract_results(results_data)
    except ValueError as e:
        print(f"ERROR extracting results: {e}")
        sys.exit(1)

    print(f"  ✓ lam_min = {results_dict['lam_min']}")
    print(f"  ✓ zeros count = {len(results_dict['zeros'])}")

    # Step 3: Precision analysis
    print("[3/5] Computing precision metrics...")
    extracted_zeros = {i + 1: z for i, z in enumerate(results_dict["zeros"])}
    metrics = compute_precision_metrics(extracted_zeros)
    print(f"  ✓ Analyzed {len(metrics)} zeros")

    # Step 4: Outlier detection
    print("[4/5] Detecting outliers...")
    outliers = identify_outliers(metrics, sigma_threshold=2.0)
    print(f"  ✓ Found {len(outliers)} outliers")

    # k=7 analysis
    k7_analysis = analyze_k7_outlier(extracted_zeros, results_dict["lam_min"])
    if k7_analysis:
        print(f"  ✓ k=7 analysis complete: |sin(γ₇L/2)| = {k7_analysis['abs_sin(gamma7 * L / 2)']:.3e}")

    # Step 5: Fetch chain log
    print("[5/5] Fetching chain log for verification...")
    chain_log = fetch_chain_log(OLYMPUS_HOST, REMOTE_CHAIN_LOG)
    chain_log_snippet = chain_log[:500] if chain_log else "[No chain log]"

    # Step 6: Format and append findings
    print(f"[6/6] Appending to {LOCAL_FINDINGS_FILE}")
    findings = format_findings(results_dict, metrics, outliers, k7_analysis, chain_log_snippet)

    try:
        with open(LOCAL_FINDINGS_FILE, "a") as f:
            f.write(findings)
        print(f"  ✓ Findings appended")
    except Exception as e:
        print(f"ERROR writing findings: {e}")
        sys.exit(1)

    # Summary
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"lam_min (sign): {results_dict['lam_min']} ({'+' if results_dict['lam_min'] > 0 else '-'})")
    print(f"Zeros extracted: {len(metrics)}")
    if metrics:
        precisions = [m["digit_precision"] for m in metrics]
        print(f"Mean digit precision: {sum(precisions) / len(precisions):.2f}")
    print(f"Outliers: {len(outliers)}")
    print(f"Findings written to: {LOCAL_FINDINGS_FILE}")
    print("=" * 70)


if __name__ == "__main__":
    main()
