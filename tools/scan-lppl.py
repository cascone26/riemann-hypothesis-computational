#!/usr/bin/env python3
"""
scan-lppl.py — Log-Periodic Power Law analysis on blob hump structure.
Tony Module 4 analog: tests whether blob peaks are log-periodically spaced.
If they are, ω = 2π / ln(r) gives a candidate Riemann zero frequency.

Sornette's LPPL: ln p(t) = A + B|tc−t|^m + C|tc−t|^m · cos(ω·ln|tc−t| + φ)
Riemann explicit formula: ψ(x) contains x^(1/2+it) = x^(1/2) · e^(it·ln x)
Same log-periodic structure; ω in LPPL ↔ Im(ρ) for Riemann zeros.

Usage:
  python3 scan-lppl.py /tmp/dense_sys_AT.log [--json out.json]
"""
import sys, re, json, math, argparse
from pathlib import Path

DEEP_THRESH = 1e-10
RIEMANN_ZEROS = [14.1347, 21.0220, 25.0109, 30.4249, 32.9351,
                 37.5862, 40.9187, 43.3271, 48.0052, 49.7738]  # Im(ρ₁..₁₀)

def parse_log(path):
    rows = []
    with open(path) as f:
        for line in f:
            m = re.match(r'\s+(\d+)\s+(-[\d.e+\-]+)\s+(-[\d.e+\-]+)\s+\d+s', line)
            if m:
                rows.append({'lsq': int(m.group(1)), 'l1': float(m.group(2))})
    return rows

def extract_blobs(rows):
    blobs, in_blob, blob_rows = [], False, []
    for r in rows:
        deep = abs(r['l1']) < DEEP_THRESH
        if not deep:
            blob_rows.append(r)
            in_blob = True
        else:
            if in_blob and len(blob_rows) >= 8:
                blobs.append(blob_rows[:])
            in_blob, blob_rows = False, []
    if in_blob and len(blob_rows) >= 8:
        blobs.append(blob_rows)
    return blobs

def find_humps(l1_vals, lsqs, min_window=3):
    """Find local maxima (humps) in blob λ₁ profile."""
    n = len(l1_vals)
    window = max(min_window, n // 15)
    maxima = []
    for i in range(window, n - window):
        if all(l1_vals[i] >= l1_vals[j] for j in range(i-window, i+window+1) if j != i):
            maxima.append(i)
    return maxima

def test_log_periodicity(hump_positions, ref_lsq):
    """
    Test if hump positions follow log-periodic spacing relative to ref_lsq.
    For LPPL: peaks at tc - (tc-t0)*r^n → diffs from ref form geometric sequence.
    """
    if len(hump_positions) < 3:
        return None
    diffs = [p - ref_lsq for p in hump_positions if p > ref_lsq]
    if len(diffs) < 3:
        return None
    ratios = []
    for i in range(len(diffs)-1):
        if diffs[i] > 0:
            ratios.append(diffs[i+1] / diffs[i])
    if not ratios:
        return None
    mean_r = sum(ratios) / len(ratios)
    std_r = (sum((r - mean_r)**2 for r in ratios) / len(ratios)) ** 0.5
    cv = std_r / mean_r if mean_r > 0 else float('inf')
    is_lp = (cv < 0.35) and (1.05 < mean_r < 5.0)
    omega = (2 * math.pi / math.log(mean_r)) if (is_lp and mean_r > 1) else None
    closest_zero = None
    if omega:
        diffs_to_zeros = [(abs(omega - z), z) for z in RIEMANN_ZEROS]
        diffs_to_zeros.sort()
        d, z = diffs_to_zeros[0]
        if d < 3.0:
            closest_zero = {'im_rho': z, 'delta': round(d, 4)}
    return {
        'n_humps': len(diffs),
        'diffs_from_start': [int(d) for d in diffs],
        'ratios': [round(r, 4) for r in ratios],
        'mean_r': round(mean_r, 4),
        'cv': round(cv, 4),
        'is_log_periodic': is_lp,
        'omega': round(omega, 4) if omega else None,
        'closest_riemann_zero': closest_zero,
    }

def analyze_blob(blob_rows, blob_idx):
    lsqs = [r['lsq'] for r in blob_rows]
    l1s  = [abs(r['l1']) for r in blob_rows]
    peak_idx = l1s.index(max(l1s))
    hump_idxs = find_humps(l1s, lsqs)
    hump_lsqs = [lsqs[i] for i in hump_idxs]
    lp = test_log_periodicity(hump_lsqs, lsqs[0])
    return {
        'blob': blob_idx + 1,
        'start': lsqs[0],
        'end': lsqs[-1],
        'span': len(blob_rows),
        'peak_l1': round(max(l1s), 4),
        'peak_lsq': lsqs[peak_idx],
        'hump_count': len(hump_idxs),
        'hump_lsqs': hump_lsqs,
        'log_periodicity': lp,
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('log', help='Scan log path')
    ap.add_argument('--json', help='JSON output path')
    args = ap.parse_args()

    rows = parse_log(args.log)
    if not rows:
        print('No data'); sys.exit(1)

    blobs = extract_blobs(rows)
    print(f'{Path(args.log).name}: {len(blobs)} blobs')
    print(f'Known Riemann zero Im(ρ): {RIEMANN_ZEROS[:5]}...')

    results = []
    lp_count = 0

    for i, blob in enumerate(blobs):
        res = analyze_blob(blob, i)
        results.append(res)
        lp = res['log_periodicity']
        flag = ''
        if lp and lp['is_log_periodic']:
            lp_count += 1
            flag = f'  *** LOG-PERIODIC r={lp["mean_r"]:.3f} ω={lp["omega"]}'
            if lp['closest_riemann_zero']:
                z = lp['closest_riemann_zero']
                flag += f'  ≈ Im(ρ)={z["im_rho"]} (Δ={z["delta"]})'
        print(f'\nBlob {i+1}: lsq {res["start"]}–{res["end"]} ({res["span"]} lsq) | '
              f'peak={res["peak_l1"]:.3e} at {res["peak_lsq"]}{flag}')
        print(f'  {res["hump_count"]} humps at {res["hump_lsqs"]}')
        if lp:
            print(f'  ratios={lp["ratios"]} mean_r={lp["mean_r"]} CV={lp["cv"]}')

    print(f'\n{lp_count}/{len(blobs)} blobs show log-periodic hump structure')
    if lp_count:
        print('→ Log-periodic blobs are consistent with the Sornette-Riemann connection:')
        print('  blob hump spacing r ↔ crash acceleration → same math as Riemann zero frequency ω')

    out = {'log': str(args.log), 'n_blobs': len(blobs), 'n_log_periodic': lp_count, 'blobs': results}
    if args.json:
        Path(args.json).write_text(json.dumps(out, indent=2))
        print(f'JSON → {args.json}')

if __name__ == '__main__':
    main()
