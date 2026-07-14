#!/usr/bin/env python3
"""
scan-phi.py — Test BC inter-arrival gaps and zone/blob spans for φ-ratio clustering.
Tony's 8/6/5 lattice: if primes carry the 5-fold aperiodic (φ) signature, φ-ratios
should appear in the gap structure of the Connes-Weil eigenvalue scan.

Tests:
  1. BC gap ratios: consecutive gaps g[i+1]/g[i] — cluster near φ, φ², 1/φ?
  2. Zone span ratios: consecutive zone spans — φ signature?
  3. Blob peak ratios: consecutive blob peaks — geometric φ progression?
  4. E-spacing test (buzsaki analog): are log-gaps spaced by ln(φ) vs ln(e) vs ln(2)?

Usage:
  python3 scan-phi.py /tmp/riemann_events.db [--json out.json]
"""
import sys, json, math, argparse, sqlite3
from pathlib import Path

PHI   = (1 + math.sqrt(5)) / 2   # 1.6180...
PHI2  = PHI ** 2                  # 2.6180...
IPHI  = 1 / PHI                   # 0.6180...
IPHI2 = 1 / PHI2                  # 0.3820...
SQRT2 = math.sqrt(2)              # silver ratio analog
SQRT3 = math.sqrt(3)              # bronze ratio analog

METALLIC = {
    '1/φ²': IPHI2,
    '1/φ':  IPHI,
    '1':    1.0,
    'φ':    PHI,
    'φ²':   PHI2,
    'e':    math.e,
    '2':    2.0,
    'π':    math.pi,
}

TOL = 0.08  # within 8% of a metallic ratio = hit

def nearest_metallic(r):
    best_name, best_dist = None, float('inf')
    for name, val in METALLIC.items():
        d = abs(r - val) / val
        if d < best_dist:
            best_dist = d
            best_name = name
    return best_name, best_dist

def ratio_histogram(ratios, label):
    hits = {name: 0 for name in METALLIC}
    misses = 0
    for r in ratios:
        name, dist = nearest_metallic(r)
        if dist < TOL:
            hits[name] += 1
        else:
            misses += 1
    total = len(ratios)
    print(f'\n{label} ({total} ratios):')
    for name, count in sorted(hits.items(), key=lambda x: -x[1]):
        if count > 0:
            pct = count / total * 100
            bar = '█' * int(pct / 2)
            print(f'  {name:6s}: {count:3d} ({pct:5.1f}%) {bar}')
    print(f'  other : {misses:3d} ({misses/total*100:5.1f}%)')
    best = max(hits, key=hits.get)
    return hits, best

def e_spacing_test(values, label):
    """
    Buzsaki analog: test whether log-spacings cluster at ln(φ), ln(e), ln(2), ln(π).
    """
    if len(values) < 4:
        return None
    log_vals = sorted([math.log(v) for v in values if v > 0])
    gaps = [log_vals[i+1] - log_vals[i] for i in range(len(log_vals)-1)]
    if not gaps:
        return None
    mean_gap = sum(gaps) / len(gaps)
    targets = {
        'ln(φ)=0.481': math.log(PHI),
        'ln(√2)=0.347': math.log(SQRT2),
        'ln(2)=0.693':  math.log(2),
        'ln(e)=1.000':  1.0,
        'ln(π)=1.145':  math.log(math.pi),
    }
    print(f'\nE-spacing test on {label} (mean log-gap={mean_gap:.4f}):')
    closest_name = min(targets, key=lambda k: abs(mean_gap - targets[k]))
    for name, val in targets.items():
        diff = abs(mean_gap - val)
        flag = ' ←' if name == closest_name else ''
        print(f'  {name}: dist={diff:.4f}{flag}')
    return closest_name, mean_gap

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('db', help='SQLite DB path')
    ap.add_argument('--json', help='JSON output path')
    args = ap.parse_args()

    conn = sqlite3.connect(args.db)

    bc_lsqs  = [r[0] for r in conn.execute('SELECT lsq FROM bc_events ORDER BY lsq').fetchall()]
    zones    = conn.execute('SELECT start_lsq, end_lsq, span FROM zones ORDER BY start_lsq').fetchall()
    blobs    = conn.execute('SELECT start_lsq, end_lsq, span, peak_l1 FROM blobs ORDER BY start_lsq').fetchall()
    conn.close()

    print(f'Dataset: {len(bc_lsqs)} BCs | {len(zones)} zones | {len(blobs)} blobs')
    print(f'φ = {PHI:.6f} | 1/φ = {IPHI:.6f} | φ² = {PHI2:.6f}')
    print(f'Tolerance: ±{TOL*100:.0f}% of each ratio')

    results = {}

    # 1. BC inter-arrival gap ratios
    if len(bc_lsqs) >= 4:
        gaps = [bc_lsqs[i+1] - bc_lsqs[i] for i in range(len(bc_lsqs)-1)]
        ratios = [gaps[i+1] / gaps[i] for i in range(len(gaps)-1) if gaps[i] > 0]
        hits, best = ratio_histogram(ratios, 'BC gap ratios g[i+1]/g[i]')
        results['bc_gaps'] = {'best': best, 'hits': hits, 'n': len(ratios)}
        e_spacing_test(gaps, 'BC gaps')

    # 2. Zone span ratios
    if len(zones) >= 4:
        spans = [z[2] for z in zones if z[2] > 0]
        ratios = [spans[i+1] / spans[i] for i in range(len(spans)-1)]
        hits, best = ratio_histogram(ratios, 'Zone span ratios')
        results['zone_spans'] = {'best': best, 'hits': hits, 'n': len(ratios)}
        e_spacing_test(spans, 'zone spans')

    # 3. Blob peak magnitude ratios
    if len(blobs) >= 4:
        peaks = [abs(b[3]) for b in blobs if b[3] and abs(b[3]) > 0]
        ratios = [peaks[i+1] / peaks[i] for i in range(len(peaks)-1)]
        hits, best = ratio_histogram(ratios, 'Blob peak magnitude ratios')
        results['blob_peaks'] = {'best': best, 'hits': hits, 'n': len(ratios)}
        e_spacing_test(peaks, 'blob peak magnitudes')

    # 4. Zone start positions: are zones themselves φ-spaced in lsq?
    if len(zones) >= 4:
        zone_starts = [z[0] for z in zones]
        gaps = [zone_starts[i+1] - zone_starts[i] for i in range(len(zone_starts)-1)]
        ratios = [gaps[i+1] / gaps[i] for i in range(len(gaps)-1) if gaps[i] > 0]
        hits, best = ratio_histogram(ratios, 'Zone position gap ratios')
        results['zone_positions'] = {'best': best, 'hits': hits, 'n': len(ratios)}

    # Summary
    print('\n=== φ SIGNATURE SUMMARY ===')
    phi_dominant = sum(1 for v in results.values() if v['best'] in ('φ', '1/φ', 'φ²', '1/φ²'))
    print(f'φ-family dominant in {phi_dominant}/{len(results)} test categories')
    for cat, v in results.items():
        phi_hit = sum(v['hits'].get(k, 0) for k in ('φ', '1/φ', 'φ²', '1/φ²'))
        total = v['n']
        print(f'  {cat}: best={v["best"]} | φ-family={phi_hit}/{total} ({phi_hit/total*100:.1f}%)')

    if args.json:
        out = {k: {kk: (vv if not isinstance(vv, dict) else {kkk: int(vvv) for kkk, vvv in vv.items()})
                   for kk, vv in v.items()} for k, v in results.items()}
        Path(args.json).write_text(json.dumps(out, indent=2))
        print(f'\nJSON → {args.json}')

if __name__ == '__main__':
    main()
