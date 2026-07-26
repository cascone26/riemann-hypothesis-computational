#!/usr/bin/env python3
"""Parse the output log from extract_ldlt_scaling_adaptive.py and extract results."""

import re
import sys
import json
from pathlib import Path

def parse_log(logfile):
    """Parse log output and extract c, lambda_min results."""
    results = []
    
    with open(logfile, 'r') as f:
        content = f.read()
    
    # Pattern: [RUN] c=XX... followed by ✓ λ_min = Y.YYYe±ZZ
    # We want to match blocks like:
    # [RUN] c=2, N=50, prec=600
    # ...
    # ✓ λ_min = 2.152210e-02  n_pos=101, n_neg=0
    
    pattern = r'\[RUN\] c=(\d+),.*?✓ λ_min = ([\d.e+-]+)'
    
    matches = re.findall(pattern, content, re.DOTALL)
    
    for c_str, lmin_str in matches:
        c = int(c_str)
        lambda_min = float(lmin_str)
        
        # Compute log values
        import math
        log10_lmin = math.log10(lambda_min) if lambda_min > 0 else None
        log_lmin = math.log(lambda_min) if lambda_min > 0 else None
        
        results.append({
            "c": c,
            "lambda_min": lambda_min,
            "log10_lambda_min": log10_lmin,
            "log_lambda_min": log_lmin
        })
    
    return results

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: parse_ldlt_log.py <logfile> [output.json]")
        sys.exit(1)
    
    logfile = sys.argv[1]
    output = sys.argv[2] if len(sys.argv) > 2 else "parsed_ldlt_results.json"
    
    results = parse_log(logfile)
    
    with open(output, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"Extracted {len(results)} results")
    print("\nc    | λ_min        | log10(λ_min)")
    print("-----|--------------|----------")
    for r in results[:35]:
        print(f"{r['c']:4d} | {r['lambda_min']:.6e} | {r['log10_lambda_min']:>11.2f}")
    
    if len(results) > 35:
        print(f"... and {len(results) - 35} more")
