#!/usr/bin/env python3
"""
scan-zero.py — Test whether BC event lsq positions correlate with Riemann zero structure.
Also collects all LPPL omega values across scans and matches to known zeros.

The Weil explicit formula: contributions of zero rho=1/2+it oscillate as x^(1/2)*cos(t*ln x).
If lsq plays the role of log(x), BC clusters should appear at spacings 2pi/t_n.

Usage:
  python3 scan-zero.py /tmp/riemann_events.db [--lppl-dir /tmp] [--json out.json]
"""
import sys, json, math, argparse, sqlite3, os, glob
from pathlib import Path

# First 100 Riemann zero imaginary parts (LMFDB verified)
ZEROS = [
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
    52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
    67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
    79.337375, 82.910381, 84.735493, 87.425275, 88.809111,
    92.491899, 94.651344, 95.870634, 98.831194, 101.317851,
    103.725538, 105.446623, 107.168611, 111.029536, 111.874659,
    114.320220, 116.226680, 118.790783, 121.370125, 122.946829,
    124.256819, 127.516683, 129.578704, 131.087688, 133.497737,
    134.756510, 138.116042, 139.736209, 141.123707, 143.111846,
    146.000982, 147.422765, 150.053521, 150.925257, 153.024693,
    156.112909, 157.597591, 158.849988, 161.188964, 163.030709,
    165.537069, 167.184439, 169.094515, 169.911976, 173.411536,
    174.754191, 176.441434, 178.377407, 179.916484, 182.207078,
    184.874467, 185.598783, 187.228922, 189.416158, 192.026656,
    193.079726, 195.265396, 196.876481, 198.015309, 201.264751,
    202.493594, 204.189671, 205.394698, 207.906258, 209.576509,
    211.690862, 213.347919, 214.547044, 216.169538, 219.067596,
    220.714918, 221.430705, 224.007000, 224.983324, 227.421444,
    229.337413, 231.250188, 231.987235, 233.693404, 236.524230,
]

def nearest_zero(omega, tol=2.0):
    best = min(ZEROS, key=lambda z: abs(omega - z))
    dist = abs(omega - best)
    return best, dist, dist < tol

def lsq_spacing_test(bc_lsqs):
    """
    Test: consecutive BC gaps — do they follow 2π/t_n spacing patterns?
    For zero t_n, expected gap in lsq = 2π/t_n * (some scale factor).
    We solve for scale factor from data, then test fit.
    """
    if len(bc_lsqs) < 10:
        return None
    gaps = [bc_lsqs[i+1] - bc_lsqs[i] for i in range(len(bc_lsqs)-1)]
    mean_gap = sum(gaps) / len(gaps)
    # If gap = 2π/t * scale, then scale = gap * t / 2π
    # For each zero, compute implied scale and variance
    results = []
    for t in ZEROS[:20]:
        expected_period = 2 * math.pi / t
        implied_scales = [g / expected_period for g in gaps]
        mean_scale = sum(implied_scales) / len(implied_scales)
        var_scale = sum((s - mean_scale)**2 for s in implied_scales) / len(implied_scales)
        cv = var_scale**0.5 / mean_scale if mean_scale > 0 else float('inf')
        results.append({'zero': t, 'mean_scale': round(mean_scale, 2), 'cv': round(cv, 3)})
    results.sort(key=lambda x: x['cv'])
    return results[:5]

def lppl_omega_survey(lppl_dir):
    """Collect all omega values from LPPL JSON files."""
    omegas = []
    for path in sorted(glob.glob(os.path.join(lppl_dir, 'dense_sys_*_lppl.json'))):
        try:
            data = json.loads(Path(path).read_text())
            label = Path(path).stem.replace('dense_sys_', '').replace('_lppl', '')
            for blob in data.get('blobs', []):
                lp = blob.get('log_periodicity')
                if lp and lp.get('is_log_periodic') and lp.get('omega'):
                    omega = lp['omega']
                    zero, dist, hit = nearest_zero(omega)
                    omegas.append({
                        'label': label,
                        'blob_start': blob['start'],
                        'omega': omega,
                        'r': lp['mean_r'],
                        'nearest_zero': zero,
                        'dist': round(dist, 4),
                        'hit': hit,
                    })
        except Exception:
            pass
    return omegas

def bc_zone_correlation(bc_lsqs, zones):
    """
    For each Riemann zero t_n, check if BC density oscillates with period 2π/t_n.
    Simple: bin BC events, compute DFT of BC density, find peaks.
    """
    if len(bc_lsqs) < 20:
        return None
    lsq_min, lsq_max = bc_lsqs[0], bc_lsqs[-1]
    span = lsq_max - lsq_min
    bin_size = max(10, span // 100)
    n_bins = span // bin_size
    bins = [0] * n_bins
    for lsq in bc_lsqs:
        b = min(int((lsq - lsq_min) / bin_size), n_bins - 1)
        bins[b] += 1
    # DFT — find dominant frequency
    freqs = []
    for k in range(1, n_bins // 2):
        re = sum(bins[n] * math.cos(2*math.pi*k*n/n_bins) for n in range(n_bins))
        im = sum(bins[n] * math.sin(2*math.pi*k*n/n_bins) for n in range(n_bins))
        power = re**2 + im**2
        freq_lsq = k / span  # cycles per lsq
        freq_t = freq_lsq * 2 * math.pi  # angular frequency
        freqs.append((power, freq_t, k))
    freqs.sort(reverse=True)
    top = freqs[:10]
    results = []
    for power, freq_t, k in top:
        zero, dist, hit = nearest_zero(freq_t)
        results.append({
            'freq_t': round(freq_t, 4),
            'power': round(power, 2),
            'nearest_zero': zero,
            'dist': round(dist, 4),
            'hit': hit,
        })
    return results

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('db', help='SQLite DB path')
    ap.add_argument('--lppl-dir', default='/tmp', help='Directory with *_lppl.json files')
    ap.add_argument('--json', help='JSON output path')
    args = ap.parse_args()

    conn = sqlite3.connect(args.db)
    bc_lsqs = [r[0] for r in conn.execute('SELECT lsq FROM bc_events ORDER BY lsq').fetchall()]
    zones = conn.execute('SELECT start_lsq, end_lsq, span FROM zones ORDER BY start_lsq').fetchall()
    conn.close()

    print(f'BC events: {len(bc_lsqs)} | lsq {bc_lsqs[0]}–{bc_lsqs[-1]}')
    print(f'First 10 Riemann zeros: {ZEROS[:10]}')

    results = {}

    # 1. LPPL omega survey
    print('\n=== LPPL ω VALUES vs RIEMANN ZEROS ===')
    omegas = lppl_omega_survey(args.lppl_dir)
    hits = [o for o in omegas if o['hit']]
    print(f'Log-periodic blobs found: {len(omegas)} | Zero hits (Δ<2): {len(hits)}')
    for o in omegas:
        flag = ' *** HIT' if o['hit'] else ''
        print(f'  Scan {o["label"]} blob@{o["blob_start"]}: ω={o["omega"]:.4f} r={o["r"]:.4f} → nearest ρ={o["nearest_zero"]} (Δ={o["dist"]:.4f}){flag}')
    results['lppl_omegas'] = omegas

    # 2. BC gap vs zero period test
    print('\n=== BC GAP SPACING vs ZERO PERIODS ===')
    spacing = lsq_spacing_test(bc_lsqs)
    if spacing:
        print('Best-fit zeros (lowest CV of gap/period ratio):')
        for r in spacing:
            print(f'  t={r["zero"]:.3f}: scale={r["mean_scale"]:.1f} lsq/period, CV={r["cv"]:.3f}')
    results['spacing_test'] = spacing

    # 3. DFT of BC density
    print('\n=== DFT OF BC EVENT DENSITY ===')
    dft = bc_zone_correlation(bc_lsqs, zones)
    if dft:
        print('Top frequency components vs Riemann zeros:')
        for r in dft[:5]:
            flag = ' *** ZERO HIT' if r['hit'] else ''
            print(f'  ω={r["freq_t"]:.4f} (power={r["power"]:.0f}) → nearest ρ={r["nearest_zero"]} (Δ={r["dist"]:.4f}){flag}')
    results['dft'] = dft

    if args.json:
        Path(args.json).write_text(json.dumps(results, indent=2))
        print(f'\nJSON → {args.json}')

if __name__ == '__main__':
    main()
