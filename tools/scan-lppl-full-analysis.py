#!/usr/bin/env python3
"""
Full LPPL spectral analysis: omega histogram, clustering test, phase coherence.
Usage: python3 scan-lppl-full-analysis.py /tmp/lppl_results_all.json

Tests:
1. Omega histogram — visual distribution vs Riemann zero positions
2. Clustering test — binomial test for enrichment near zeros
3. Phase coherence — Rayleigh test for global phase lock

The theoretical connection:
  If λ₁(Λ) ~ (tc-Λ)^(½+i·Im(ρ₁)) near blob exit tc,
  then LPPL fit recovers ω = Im(ρ₁) = 14.1347.
  Phase coherence (φ - Im(ρ₁)·ln(tc) = const) would confirm
  that oscillations are driven by a SINGLE global spectral frequency.
"""
import json, math, sys, os
from math import comb, log, exp, lgamma
from collections import defaultdict

ZEROS = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178,
         40.918719, 43.327073, 48.005151, 49.773832]
ZLBLS = ['rho1','rho2','rho3','rho4','rho5','rho6','rho7','rho8','rho9','rho10']

def tval_for_df(df):
    t = {1:12.71,2:4.303,3:3.182,4:2.776,5:2.571,6:2.447,7:2.365,8:2.306,9:2.262,10:2.228,11:2.201,12:2.179}
    return t.get(df, 2.0)

def stats(ws, zi):
    n = len(ws)
    if n == 0: return None
    mean = sum(ws)/n
    if n > 1:
        sd = math.sqrt(sum((x-mean)**2 for x in ws)/(n-1))
        se = sd/math.sqrt(n)
        t = (mean-zi)/se
        tv = tval_for_df(n-1)
        return {'n':n,'mean':mean,'delta':mean-zi,'se':se,'t':t,'lo':mean-tv*se,'hi':mean+tv*se}
    return {'n':1,'mean':mean,'delta':mean-zi,'se':None,'t':None,'lo':None,'hi':None}

def binom_pvalue_ge(k, n, p):
    if k <= 0: return 1.0
    prob = 0.0
    for i in range(k, n+1):
        log_term = lgamma(n+1)-lgamma(i+1)-lgamma(n-i+1) + i*log(p) + (n-i)*log(max(1-p,1e-300))
        prob += exp(log_term)
    return min(1.0, prob)

jsonfile = sys.argv[1] if len(sys.argv) > 1 else '/tmp/lppl_results_all.json'
data = json.load(open(jsonfile))
hits = [h for h in data if h.get('hit')]
fitted = [h for h in data if 'omega' in h and h['omega'] is not None]

print(f"=== LPPL Full Analysis: {os.path.basename(jsonfile)} ===")
print(f"Total blobs: {len(data)}, fitted: {len(fitted)}, hits (delta<2): {len(hits)}")

# === Section 1: Zero cluster statistics ===
print("\n--- 1. Zero cluster statistics ---")
by_zero = defaultdict(list)
for h in hits:
    by_zero[h['nearest_zero']].append(h['omega'])

for zi, lbl in zip(ZEROS, ZLBLS):
    ws = by_zero.get(zi, [])
    if not ws: continue
    s = stats(ws, zi)
    r_theory = math.exp(2*math.pi/zi)
    implied_rs = [math.exp(2*math.pi/w) for w in ws]
    mean_r = sum(implied_rs)/len(implied_rs)
    if s['se'] is not None:
        print(f"  {lbl}: n={s['n']} mean_w={s['mean']:.4f} delta={s['delta']:+.4f} SE={s['se']:.4f} t={s['t']:.3f} 95%CI=[{s['lo']:.3f},{s['hi']:.3f}]")
        print(f"         implied_r mean={mean_r:.4f} theory={r_theory:.4f} err={((mean_r-r_theory)/r_theory)*100:+.2f}%")
    else:
        print(f"  {lbl}: n=1 w={ws[0]:.4f} delta={s['delta']:+.4f} implied_r={math.exp(2*math.pi/ws[0]):.4f} theory={r_theory:.4f}")

# === Section 2: Weil amplitude chi2 ===
print("\n--- 2. Weil amplitude chi2 (rho1:rho2:rho3) ---")
counts = [len(by_zero.get(zi,[])) for zi in ZEROS[:3]]
pred_ratio = [1/z for z in ZEROS[:3]]
total_pred = sum(pred_ratio)
total_obs = sum(counts)
expected = [r/total_pred*total_obs for r in pred_ratio]
chi2 = sum((o-e)**2/e for o,e in zip(counts,expected) if e>0)
print(f"  Observed:  {' '.join(f'{lbl}={c}' for lbl,c in zip(ZLBLS[:3],counts))}")
print(f"  Expected:  {' '.join(f'{lbl}={e:.2f}' for lbl,e in zip(ZLBLS[:3],expected))}")
print(f"  chi2={chi2:.4f} (2 dof)")

# === Section 3: Omega histogram ===
print("\n--- 3. Omega histogram ---")
omegas_all = [h['omega'] for h in fitted]
BIN = 0.5
lo, hi = 10.0, 47.0
bins = {}
for o in omegas_all:
    b = int((o-lo)/BIN)
    bins[b] = bins.get(b,0)+1
zero_bins = {int((z-lo)/BIN): lbl for z,lbl in zip(ZEROS,ZLBLS) if lo<=z<=hi}
n_bins = int((hi-lo)/BIN)
max_count = max(bins.values()) if bins else 1
scale = 30/max_count
for b in range(n_bins):
    bstart = lo + b*BIN
    cnt = bins.get(b, 0)
    bar = '#' * int(cnt*scale)
    zero_mark = f" <-- {zero_bins[b]}" if b in zero_bins else ""
    if cnt > 0 or b in zero_bins:
        print(f"  {bstart:5.1f}: {bar:<30} ({cnt:2d}){zero_mark}")

# === Section 4: Clustering test ===
print("\n--- 4. Clustering enrichment test ---")
DELTA = 1.0
near_zero_count = 0
for zi, lbl in zip(ZEROS, ZLBLS):
    n_near = sum(1 for o in omegas_all if abs(o-zi) < DELTA)
    if n_near > 0:
        print(f"  {lbl}: {n_near} blobs within delta<{DELTA}")
    near_zero_count += n_near

zeros_in_range = [z for z in ZEROS if lo<=z<=hi]
total_zero_width = len(zeros_in_range) * 2*DELTA
total_range = hi - lo
expected_frac = total_zero_width / total_range
expected_near = len(omegas_all) * expected_frac
print(f"  Total near zeros: {near_zero_count}/{len(omegas_all)} (expected {expected_near:.1f} = {expected_frac*100:.1f}% of range)")
print(f"  Enrichment: {near_zero_count/expected_near:.2f}x")
p_val = binom_pvalue_ge(near_zero_count, len(omegas_all), expected_frac)
print(f"  Binomial p-value: {p_val:.6f}")

# === Section 5: Phase coherence ===
print("\n--- 5. Phase coherence (rho1) ---")
rho1_hits = [h for h in hits if h['nearest_zero'] == 14.134725 and 'phi' in h]
print(f"  rho1 hits with phase: {len(rho1_hits)}")
if len(rho1_hits) >= 2:
    global_phis = []
    for h in rho1_hits:
        phi = h['phi']
        tc = h.get('tc', h['blob_end']+1)
        gp = (phi - 14.134725 * math.log(tc)) % (2*math.pi)
        global_phis.append(gp)
        lbl = h.get('label','?')
        print(f"  {h['blob_start']:>6}@{lbl}: phi={phi:.3f} tc={tc} global_phi={gp:.3f}")
    c = sum(math.cos(a) for a in global_phis)/len(global_phis)
    s = sum(math.sin(a) for a in global_phis)/len(global_phis)
    R = math.sqrt(c**2+s**2)
    mean_angle = math.atan2(s,c) % (2*math.pi)
    n = len(global_phis)
    p_rayleigh = math.exp(-n * R**2) if n > 0 else 1.0
    print(f"\n  Rayleigh R = {R:.4f} (0=random, 1=coherent)")
    print(f"  Mean global phase = {mean_angle:.4f} rad = {math.degrees(mean_angle):.1f} deg")
    print(f"  Rayleigh p-value ≈ {p_rayleigh:.4f}")
