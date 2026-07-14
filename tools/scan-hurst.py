#!/usr/bin/env python3
"""
scan-hurst.py — Hurst exponent on λ₁(lsq) eigenvalue sequence.
Tony Module 2/5 analog: applies fractal / multifractal analysis to Riemann spectral data.

Usage:
  python3 scan-hurst.py /tmp/dense_sys_AT.log [--json out.json]
  python3 scan-hurst.py /tmp/dense_sys_AT.log /tmp/dense_sys_AU.log  # concatenate
"""
import sys, re, json, math, argparse
from pathlib import Path

def parse_log(path):
    rows = []
    with open(path) as f:
        for line in f:
            m = re.match(r'\s+(\d+)\s+(-[\d.e+\-]+)\s+(-[\d.e+\-]+)\s+\d+s', line)
            if m:
                rows.append((int(m.group(1)), float(m.group(2))))
    return rows

def rs_analysis(series):
    """R/S (rescaled range) analysis for Hurst exponent."""
    n = len(series)
    sizes, rs_vals = [], []
    size = 16
    while size <= n // 2:
        chunks = []
        for start in range(0, n - size, size):
            chunk = series[start:start+size]
            mean = sum(chunk) / len(chunk)
            devs = [x - mean for x in chunk]
            cumdev, c = [], 0
            for d in devs:
                c += d
                cumdev.append(c)
            R = max(cumdev) - min(cumdev)
            S = (sum((x - mean)**2 for x in chunk) / len(chunk)) ** 0.5
            if S > 0:
                chunks.append(R / S)
        if chunks:
            sizes.append(size)
            rs_vals.append(sum(chunks) / len(chunks))
        size = min(size * 2, size + 64) if size < 128 else size * 2
    if len(sizes) < 3:
        return None
    log_s = [math.log(s) for s in sizes]
    log_r = [math.log(r) for r in rs_vals]
    n_pts = len(log_s)
    sx = sum(log_s); sy = sum(log_r)
    sxy = sum(x*y for x,y in zip(log_s, log_r))
    sxx = sum(x**2 for x in log_s)
    return (n_pts*sxy - sx*sy) / (n_pts*sxx - sx**2)

def dfa(series):
    """Detrended Fluctuation Analysis."""
    n = len(series)
    mean = sum(series) / n
    profile, c = [], 0
    for x in series:
        c += x - mean
        profile.append(c)
    sizes, F_vals = [], []
    size = 16
    while size <= n // 4:
        flucts = []
        for start in range(0, n - size, size):
            chunk = profile[start:start+size]
            xm = (size - 1) / 2
            ym = sum(chunk) / size
            num = sum((i - xm)*(chunk[i] - ym) for i in range(size))
            den = sum((i - xm)**2 for i in range(size))
            if den > 0:
                slope = num / den
                intercept = ym - slope * xm
                res = [chunk[i] - (slope*i + intercept) for i in range(size)]
                flucts.append(sum(r**2 for r in res) / size)
        if flucts:
            F = (sum(flucts) / len(flucts)) ** 0.5
            if F > 0:
                sizes.append(size)
                F_vals.append(F)
        size = int(size * 1.5) if size < 128 else size + 128
    if len(sizes) < 3:
        return None
    log_s = [math.log(s) for s in sizes]
    log_F = [math.log(f) for f in F_vals]
    n_pts = len(log_s)
    sx = sum(log_s); sy = sum(log_F)
    sxy = sum(x*y for x,y in zip(log_s, log_F))
    sxx = sum(x**2 for x in log_s)
    return (n_pts*sxy - sx*sy) / (n_pts*sxx - sx**2)

def interpret_H(H):
    if H < 0.4:
        return 'anti-persistent: eigenvalue mean-reverts (zones actively balance blobs)'
    elif H < 0.6:
        return 'random walk: eigenvalue increments are memoryless'
    elif H < 0.8:
        return 'persistent: eigenvalue magnitude has memory (blobs cluster)'
    else:
        return 'strongly persistent: long-range correlations in spectral minimum'

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('logs', nargs='+', help='Scan log path(s)')
    ap.add_argument('--json', help='JSON output path')
    args = ap.parse_args()

    all_rows = []
    for log_path in args.logs:
        all_rows.extend(parse_log(log_path))
    all_rows.sort(key=lambda r: r[0])

    if not all_rows:
        print('No data found'); sys.exit(1)

    series_raw = [abs(r[1]) for r in all_rows]
    series = [math.log10(x) for x in series_raw if x > 0]

    print(f'λ₁ sequence: {len(series)} points, lsq {all_rows[0][0]}–{all_rows[-1][0]}')

    H_rs  = rs_analysis(series)
    H_dfa = dfa(series)
    H_best = H_dfa or H_rs

    print(f'\nHurst (R/S):  {H_rs:.4f}' if H_rs else '\nHurst (R/S):  insufficient data')
    print(f'Hurst (DFA):  {H_dfa:.4f}' if H_dfa else 'Hurst (DFA):  insufficient data')
    if H_best:
        print(f'\nInterpretation: {interpret_H(H_best)}')
        print(f'Tony Module 5 analog: λ₁(lsq) fractal structure {"supports mean-reversion (anti-persistent)" if H_best < 0.5 else "shows long-range memory (persistent)"}')

    result = {
        'logs': [str(p) for p in args.logs],
        'n_points': len(series),
        'lsq_range': [all_rows[0][0], all_rows[-1][0]],
        'H_rs':  round(H_rs, 4) if H_rs else None,
        'H_dfa': round(H_dfa, 4) if H_dfa else None,
        'H_best': round(H_best, 4) if H_best else None,
        'interpretation': interpret_H(H_best) if H_best else None,
    }

    if args.json:
        Path(args.json).write_text(json.dumps(result, indent=2))
        print(f'\nJSON → {args.json}')

if __name__ == '__main__':
    main()
