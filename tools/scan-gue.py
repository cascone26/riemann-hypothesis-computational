#!/usr/bin/env python3
"""
scan-gue.py — GUE spacing statistics on BC event positions.
Tests whether BC inter-arrival gaps follow GUE Wigner-Dyson vs Poisson distribution.
Montgomery-Odlyzko law: Riemann zero spacings ~ GUE eigenvalue spacings.

Usage:
  python3 scan-gue.py /tmp/riemann_events.db
  python3 scan-gue.py /tmp/dense_sys_AT_report.json [more.json ...]
"""
import sys, json, math, argparse, sqlite3
from pathlib import Path

def gue_wigner(s):
    """GUE Wigner surmise: p(s) = (32/π²) s² exp(-4s²/π)"""
    return (32 / math.pi**2) * s**2 * math.exp(-4 * s**2 / math.pi)

def goe_wigner(s):
    """GOE Wigner surmise: p(s) = (π/2) s exp(-πs²/4)"""
    return (math.pi / 2) * s * math.exp(-math.pi * s**2 / 4)

def poisson_pdf(s):
    """Poisson (uncorrelated): p(s) = exp(-s)"""
    return math.exp(-s)

def integrate_pdf(pdf, a, b, steps=500):
    """Trapezoid integration of pdf from a to b."""
    dt = (b - a) / steps
    total = 0.0
    for i in range(steps):
        t = a + i * dt
        total += (pdf(t) + pdf(t + dt)) / 2 * dt
    return total

def chi2_test(normalized_spacings, pdf, n_bins=15):
    """Chi-squared goodness-of-fit against theoretical pdf."""
    if len(normalized_spacings) < 15:
        return None
    s_max = min(max(normalized_spacings), 4.0)
    bin_w = s_max / n_bins
    observed = [0] * n_bins
    for s in normalized_spacings:
        if s < s_max:
            b = min(int(s / bin_w), n_bins - 1)
            observed[b] += 1
    n = len(normalized_spacings)
    chi2 = 0.0
    for i in range(n_bins):
        expected = integrate_pdf(pdf, i * bin_w, (i+1) * bin_w) * n
        if expected > 0.5:
            chi2 += (observed[i] - expected)**2 / expected
    return round(chi2, 3)

def load_bc_from_db(db_path):
    conn = sqlite3.connect(db_path)
    rows = conn.execute('SELECT lsq FROM bc_events ORDER BY lsq').fetchall()
    conn.close()
    return [r[0] for r in rows]

def load_bc_from_jsons(paths):
    lsqs = []
    for p in paths:
        data = json.loads(Path(p).read_text())
        for e in data.get('bc_events', []):
            lsqs.append(e['lsq'])
    return sorted(lsqs)

def analyze(lsqs):
    if len(lsqs) < 5:
        return None
    gaps = [lsqs[i+1] - lsqs[i] for i in range(len(lsqs)-1)]
    mean_gap = sum(gaps) / len(gaps)
    norm = [g / mean_gap for g in gaps]

    var = sum((s - 1.0)**2 for s in norm) / len(norm)  # mean of normalized = 1
    # level repulsion: fraction of gaps < 0.5 * mean (GUE avoids small spacings)
    repulsion = sum(1 for s in norm if s < 0.5) / len(norm)

    chi2_gue = chi2_test(norm, gue_wigner)
    chi2_goe = chi2_test(norm, goe_wigner)
    chi2_poi = chi2_test(norm, poisson_pdf)

    fits = [(name, val) for name, val in
            [('GUE (Wigner-Dyson)', chi2_gue),
             ('GOE (Wigner-Dyson)', chi2_goe),
             ('Poisson (independent)', chi2_poi)]
            if val is not None]
    fits.sort(key=lambda x: x[1])
    best = fits[0][0] if fits else 'unknown'

    return {
        'n_bc_events': len(lsqs),
        'n_gaps': len(gaps),
        'lsq_range': [lsqs[0], lsqs[-1]],
        'mean_gap_lsq': round(mean_gap, 2),
        'normalized_variance': round(var, 4),
        'expected_gue_variance': 0.286,
        'expected_poisson_variance': 1.0,
        'fraction_below_half_mean': round(repulsion, 4),
        'chi2_gue': chi2_gue,
        'chi2_goe': chi2_goe,
        'chi2_poisson': chi2_poi,
        'best_fit': best,
        'gaps': gaps,
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('inputs', nargs='+', help='SQLite DB or JSON report path(s)')
    ap.add_argument('--json', help='JSON output path')
    args = ap.parse_args()

    if args.inputs[0].endswith('.db'):
        lsqs = load_bc_from_db(args.inputs[0])
        src = args.inputs[0]
    else:
        lsqs = load_bc_from_jsons(args.inputs)
        src = f'{len(args.inputs)} JSON file(s)'

    print(f'Source: {src}')
    print(f'BC events: {len(lsqs)} | lsq {lsqs[0]}–{lsqs[-1]}')

    res = analyze(lsqs)
    if not res:
        print('Need ≥ 5 BC events for spacing analysis'); sys.exit(1)

    print(f'\nBC inter-arrival gap statistics ({res["n_gaps"]} gaps):')
    print(f'  Mean gap: {res["mean_gap_lsq"]:.1f} lsq')
    print(f'  Normalized variance: {res["normalized_variance"]:.4f}')
    print(f'    (GUE expected ~0.286 | Poisson expected 1.000)')
    print(f'  Fraction below ½·mean: {res["fraction_below_half_mean"]:.4f}')
    print(f'    (GUE level repulsion suppresses small spacings)')

    print(f'\nChi² goodness-of-fit (lower = better fit):')
    print(f'  GUE Wigner-Dyson:    χ² = {res["chi2_gue"]}')
    print(f'  GOE Wigner-Dyson:    χ² = {res["chi2_goe"]}')
    print(f'  Poisson (random):    χ² = {res["chi2_poisson"]}')
    print(f'\n→ Best fit: {res["best_fit"]}')

    if res['normalized_variance'] < 0.5:
        print('\nInterpretation: LEVEL REPULSION detected in BC spacings.')
        print('BC events avoid clustering — consistent with GUE random matrix statistics.')
        print('Aligns with Montgomery-Odlyzko: Riemann zero spacings ↔ GUE eigenvalue spacings.')
    elif res['normalized_variance'] < 0.8:
        print('\nInterpretation: Mild level repulsion — between GUE and Poisson.')
    else:
        print('\nInterpretation: Approximately Poisson — BC events are roughly independent.')

    if args.json:
        out = {k: v for k, v in res.items() if k != 'gaps'}
        out['gaps_sample'] = res['gaps'][:20]
        Path(args.json).write_text(json.dumps(out, indent=2))
        print(f'\nJSON → {args.json}')

if __name__ == '__main__':
    main()
