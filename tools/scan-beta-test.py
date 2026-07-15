#!/usr/bin/env python3
"""
Direct complex-scaling exponent test for RH.

Fits |λ₁(lsq)| = C * (tc-lsq)^β * |cos(Im(ρ₁)*ln(tc-lsq) + φ)|

for each blob that hits ρ₁ (LPPL ω within delta<2 of Im(ρ₁)=14.1347).

RH prediction: β = ½ = 0.5 (real part of zeros on critical line).
If zeros had real part σ ≠ ½, we'd see β ≠ 0.5.

This is a DIRECT test of RH via the complex scaling exponent, independent
of the LPPL fit quality or the m-exponent.

Usage: python3 scan-beta-test.py /tmp/lppl_results_all.json [log_dir]
The log_dir is needed to load the actual λ₁ data for each blob.
Default log_dir: /tmp
"""
import json, math, sys, re, glob
from pathlib import Path

RHO1 = 14.134725
RH_PREDICTION = 0.5

def load_log(log_path):
    """Load lsq -> λ₁ from a dense_sys log file."""
    rows = []
    with open(log_path) as f:
        for line in f:
            m = re.match(r'\s*(\d+)\s+(-[\d.e+\-]+)\s+', line)
            if m:
                rows.append((int(m.group(1)), float(m.group(2))))
    return rows

def fit_beta(blob_data, tc, omega=RHO1):
    """
    Fit |λ₁| = C * (tc-lsq)^β * |cos(omega*ln(tc-lsq) + φ)|
    using grid search over β and φ, linear solve for log(C).

    Returns (beta, phi, C, residual)
    """
    lsqs = [r[0] for r in blob_data]
    l1s = [abs(r[1]) for r in blob_data]

    # Filter: remove points where |λ₁| is essentially zero (would give log(-inf))
    valid = [(lsq, l1) for lsq, l1 in zip(lsqs, l1s)
             if l1 > 1e-30 and lsq < tc]
    if len(valid) < 10:
        return None

    v_lsqs, v_l1s = zip(*valid)
    log_l1s = [math.log(l1) for l1 in v_l1s]

    best_res = float('inf')
    best = None

    beta_vals = [i/20 for i in range(1, 20)]  # 0.05 to 0.95 in 0.05 steps
    phi_vals = [2*math.pi*j/32 for j in range(32)]

    for beta in beta_vals:
        for phi in phi_vals:
            # log|λ₁| = log(C) + β*log(tc-lsq) + log|cos(ω*log(tc-lsq)+φ)|
            ys = []
            xs = []
            for lsq, ll in zip(v_lsqs, log_l1s):
                u = tc - lsq
                if u <= 0: continue
                cos_term = math.cos(omega * math.log(u) + phi)
                if abs(cos_term) < 0.01:  # skip near-zeros of cosine
                    continue
                log_model_no_C = beta * math.log(u) + math.log(abs(cos_term))
                ys.append(ll)
                xs.append(log_model_no_C)

            if len(xs) < 5: continue

            # Solve for log(C): log|λ₁| = log(C) + log_model_no_C
            # log(C) = mean(log|λ₁| - log_model_no_C)
            logC = sum(y - x for y, x in zip(ys, xs)) / len(xs)

            # Residual
            res = sum((y - logC - x)**2 for y, x in zip(ys, xs)) / len(xs)

            if res < best_res:
                best_res = res
                best = (beta, phi, math.exp(logC), res)

    return best

def main():
    jsonfile = sys.argv[1] if len(sys.argv) > 1 else '/tmp/lppl_results_all.json'
    log_dir = sys.argv[2] if len(sys.argv) > 2 else '/tmp'

    data = json.load(open(jsonfile))
    rho1_hits = [h for h in data if h.get('hit') and h.get('nearest_zero') == RHO1]

    print(f"=== Direct β (complex scaling) test for RH ===")
    print(f"RH prediction: β = {RH_PREDICTION}")
    print(f"rho1 hits to test: {len(rho1_hits)}")
    print()

    betas = []

    for hit in rho1_hits:
        label = hit.get('label', '?')
        blob_start = hit['blob_start']
        blob_end = hit['blob_end']
        tc = blob_end + 1

        # Find the log file
        log_path = Path(log_dir) / f"dense_sys_{label}.log"
        if not log_path.exists():
            print(f"  {blob_start}@{label}: LOG NOT FOUND at {log_path}")
            continue

        # Load blob data
        all_rows = load_log(str(log_path))
        blob_rows = [(lsq, l1) for lsq, l1 in all_rows
                     if blob_start <= lsq <= blob_end]

        if len(blob_rows) < 10:
            print(f"  {blob_start}@{label}: TOO FEW POINTS ({len(blob_rows)})")
            continue

        result = fit_beta(blob_rows, tc, omega=RHO1)
        if result is None:
            print(f"  {blob_start}@{label}: FIT FAILED")
            continue

        beta, phi, C, res = result
        delta_from_rh = beta - RH_PREDICTION
        betas.append(beta)

        print(f"  {blob_start:>6}@{label}: β={beta:.3f} (Δ_RH={delta_from_rh:+.3f}) φ={phi:.3f} res={res:.4f} [LPPL_w={hit['omega']:.2f}]")

    if betas:
        n = len(betas)
        mean_beta = sum(betas)/n
        delta = mean_beta - RH_PREDICTION
        if n > 1:
            sd = math.sqrt(sum((b-mean_beta)**2 for b in betas)/(n-1))
            se = sd/math.sqrt(n)
            t = delta/se
            tv = {1:12.71,2:4.303,3:3.182,4:2.776,5:2.571,6:2.447,7:2.365,8:2.306}.get(n-1,2.0)
            ci = tv*se
            print(f"\n  RESULT: n={n} mean_β={mean_beta:.4f} delta={delta:+.4f} SE={se:.4f} t={t:.3f}")
            print(f"  95%CI=[{mean_beta-ci:.4f},{mean_beta+ci:.4f}]")
            print(f"  RH prediction β=0.5 {'INSIDE CI' if mean_beta-ci <= 0.5 <= mean_beta+ci else 'OUTSIDE CI'}")
        else:
            print(f"\n  RESULT: n=1 β={betas[0]:.4f} delta={betas[0]-RH_PREDICTION:+.4f}")

if __name__ == '__main__':
    main()
