"""
Fourier analysis of sigma_c(n mod T) periodicity.
Tests whether sigma_c is a well-behaved function of the phase c = n mod T.

Key questions:
1. Is sigma_c(c) bounded strictly below delta = 0.06792?
2. Does sigma_c(c) show consistent periodic structure across different k = n//T?
3. What is the Fourier decomposition? Dominant frequency, amplitude?
4. Is there a drift with k (i.e., does sigma_c change as n grows)?
"""
import csv
import numpy as np
import json

# Load data
data = []
with open('/Users/scones/projects/riemann/results/dense_sc_data.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        data.append({
            'n': int(row['n']),
            'sigma_c': float(row['sigma_c']),
            'c_mod_T': float(row['c_mod_T']),
            'k': int(float(row['k'])),
            'label': row['label'],
            'certifiable': int(row['certifiable'])
        })

delta = 0.067921
T = 88.9

print(f"Total data points: {len(data)}")
print(f"Certifiable: {sum(r['certifiable'] for r in data)}")

# Separate certifiable vs barrier
cert = [r for r in data if r['certifiable']]
barrier = [r for r in data if not r['certifiable']]

print(f"\n=== Certifiable points ===")
print(f"Count: {len(cert)}")
sc_vals = [r['sigma_c'] for r in cert]
print(f"sigma_c range: [{min(sc_vals):.6f}, {max(sc_vals):.6f}]")
print(f"delta = {delta:.6f}")
print(f"max/delta ratio: {max(sc_vals)/delta:.4f}")
print(f"Mean sigma_c: {np.mean(sc_vals):.6f}")
print(f"Std sigma_c: {np.std(sc_vals):.6f}")

# Phase analysis: c = n mod T
c_vals = np.array([r['c_mod_T'] for r in cert])
sc_arr = np.array(sc_vals)
n_arr = np.array([r['n'] for r in cert])
k_arr = np.array([r['k'] for r in cert])

print(f"\n=== Phase c = n mod T analysis ===")
print(f"c range: [{c_vals.min():.2f}, {c_vals.max():.2f}]  (T={T})")

# Bin sigma_c by phase c
n_bins = 20
bins = np.linspace(0, T, n_bins+1)
bin_centers = 0.5*(bins[:-1]+bins[1:])
bin_means = []
bin_maxes = []
bin_counts = []
for i in range(n_bins):
    mask = (c_vals >= bins[i]) & (c_vals < bins[i+1])
    if mask.sum() > 0:
        bin_means.append(np.mean(sc_arr[mask]))
        bin_maxes.append(np.max(sc_arr[mask]))
        bin_counts.append(mask.sum())
    else:
        bin_means.append(None)
        bin_maxes.append(None)
        bin_counts.append(0)

print(f"\nPhase-binned sigma_c (bins={n_bins}, T={T}):")
print(f"{'c_center':>10}  {'count':>6}  {'mean_sc':>10}  {'max_sc':>10}  {'max/delta':>10}")
for i, (bc, bm, bmx, bcnt) in enumerate(zip(bin_centers, bin_means, bin_maxes, bin_counts)):
    if bm is not None:
        print(f"{bc:10.2f}  {bcnt:6d}  {bm:10.6f}  {bmx:10.6f}  {bmx/delta:10.4f}")

# Fourier analysis on sorted phase
print(f"\n=== Fourier decomposition of sigma_c(c) ===")
# Sort by phase
sort_idx = np.argsort(c_vals)
c_sorted = c_vals[sort_idx]
sc_sorted = sc_arr[sort_idx]

# DFT: treat c as time variable in [0, T), compute FFT
# Interpolate to uniform grid for clean FFT
c_uniform = np.linspace(0, T, 256, endpoint=False)
sc_interp = np.interp(c_uniform, c_sorted, sc_sorted, period=T)

fft_vals = np.fft.rfft(sc_interp)
freqs = np.fft.rfftfreq(256, d=T/256)  # cycles per unit of c
amplitudes = np.abs(fft_vals) / 128  # normalize

print(f"DC (mean): {amplitudes[0]:.6f}")
print(f"Top 10 non-DC Fourier components:")
ndc_amps = amplitudes[1:]
ndc_freqs = freqs[1:]
top_idx = np.argsort(ndc_amps)[::-1][:10]
for idx in top_idx:
    freq = ndc_freqs[idx]
    amp = ndc_amps[idx]
    period_c = 1/freq if freq > 0 else float('inf')
    print(f"  freq={freq:.4f} (period={period_c:.2f} in c), amplitude={amp:.6f}, rel={amp/amplitudes[0]:.4f}")

# K-dependence analysis: does sigma_c drift with k?
print(f"\n=== K-drift analysis (does sigma_c change as n grows?) ===")
k_vals_cert = k_arr

# Group by k-quintile
k_sorted_idx = np.argsort(k_vals_cert)
quintile_size = len(k_vals_cert) // 5
for q in range(5):
    start = q * quintile_size
    end = (q+1)*quintile_size if q < 4 else len(k_vals_cert)
    idx_q = k_sorted_idx[start:end]
    k_q = k_vals_cert[idx_q]
    sc_q = sc_arr[idx_q]
    n_q = n_arr[idx_q]
    print(f"  Q{q+1}: n≈[{n_q.min():.0f}, {n_q.max():.0f}], "
          f"k≈[{k_q.min():.0f}, {k_q.max():.0f}], "
          f"sigma_c: mean={np.mean(sc_q):.5f}, max={np.max(sc_q):.5f}")

# Check if max sigma_c per quintile is trending up or down
q_maxes = []
for q in range(5):
    start = q * quintile_size
    end = (q+1)*quintile_size if q < 4 else len(k_vals_cert)
    idx_q = k_sorted_idx[start:end]
    q_maxes.append(np.max(sc_arr[idx_q]))
trend = np.polyfit(range(5), q_maxes, 1)
print(f"\nMax sigma_c per quintile: {[f'{x:.5f}' for x in q_maxes]}")
print(f"Linear trend slope: {trend[0]:.6f} (positive = growing, negative = shrinking)")
print(f"Trend interpretation: {'GROWING toward delta' if trend[0] > 0 else 'STABLE or shrinking'}")

# Cross-run comparison: same c_mod_T across different regions
print(f"\n=== Cross-label consistency ===")
labels = list(set(r['label'] for r in cert))
for lab in sorted(labels):
    pts = [r for r in cert if r['label'] == lab]
    sc_l = [r['sigma_c'] for r in pts]
    c_l = [r['c_mod_T'] for r in pts]
    n_l = [r['n'] for r in pts]
    print(f"  {lab}: {len(pts)} pts, n∈[{min(n_l):.0f},{max(n_l):.0f}], "
          f"sigma_c: mean={np.mean(sc_l):.5f}, max={np.max(sc_l):.5f}, "
          f"c_mod_T: [{min(c_l):.1f},{max(c_l):.1f}]")

# Save results for Opus reasoning
results = {
    'total_cert': len(cert),
    'delta': delta,
    'T': T,
    'sc_min': float(min(sc_vals)),
    'sc_max': float(max(sc_vals)),
    'sc_mean': float(np.mean(sc_vals)),
    'sc_std': float(np.std(sc_vals)),
    'max_over_delta': float(max(sc_vals)/delta),
    'fft_dc': float(amplitudes[0]),
    'fft_top_freq': float(ndc_freqs[top_idx[0]]),
    'fft_top_amp': float(ndc_amps[top_idx[0]]),
    'fft_top_rel': float(ndc_amps[top_idx[0]]/amplitudes[0]),
    'k_drift_slope': float(trend[0]),
    'q_maxes': [float(x) for x in q_maxes],
    'phase_bin_maxes': [float(x) if x else None for x in bin_maxes],
    'phase_bin_centers': [float(x) for x in bin_centers],
}
with open('/Users/scones/projects/riemann/results/fourier_sc_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to results/fourier_sc_results.json")
print(f"\n=== KEY FINDING ===")
print(f"max sigma_c = {max(sc_vals):.6f} vs delta = {delta:.6f}")
print(f"Safety margin: {delta - max(sc_vals):.6f} ({(1 - max(sc_vals)/delta)*100:.1f}%)")
if max(sc_vals) < delta:
    print("ALL CERTIFIABLE POINTS HAVE sigma_c < delta. ✓")
else:
    print("WARNING: some certifiable points have sigma_c >= delta!")
