#!/usr/bin/env python3
"""CONCENTRATION OPERATOR K = B_W · 1_{σ<0} · B_W (§888). The largest eigenvalue λ_max(K) is the
maximum fraction of |φ|² mass a bandwidth-4W Fejér function can concentrate on the symbol's negative
set {σ<0}. The razor ⟺ margin = floor − λ_max(K)·|min σ| ⪰ 0 — can the worst-case negative mass,
at its worst, beat the positive mean? Track λ_max(K) and the gap across heights. Does λ_max sit on the
Landau-Widom plunge curve (near 1/2 at the critical time-bandwidth product, and rising/falling in
robust region vs corner)? Build K in the Fejér basis, compute top eigenvalues and the positivity
certificate margin. Files: operators/concentration_operator.py"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
import json
import os
np.set_printoptions(precision=4, suppress=True, linewidth=160)

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam

print("DEBUG: Starting...", flush=True)
NMAX = 3000; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1); LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)
print(f"DEBUG: Precomputed von Mangoldt up to {NMAX}", flush=True)

def run_concentration(a, b, fourW, N, dr=0.01):
    """Build the concentration operator K = P.T @ D_σ @ P in coefficient space (Fejér basis).
    K has eigenvalues λ with λ_max = max fraction of |φ|² that can sit on {σ<0}.

    Uses memory-efficient approach: compute K iteratively without storing full P matrix.
    """
    W = fourW / 4
    # Reduced grid for memory efficiency
    RGRID = np.arange(-(b+50), (b+50), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    CEN = np.linspace(a, b, N)

    # Fejér basis: fej(x) = (sin(Wx)/(Wx))^2
    def fej(x):
        y = W * x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz]) / y[nz])**2
        return out

    # Szegő symbol σ(r) computed first
    floor_r = OMEGA / (2*np.pi)
    mask = LOGN < fourW
    ln = LOGN[mask]; lw = LAMW[mask]
    cosu = np.cos(np.outer(ln, RGRID))
    P_prime = (1/np.pi) * (lw[:, None] * cosu).sum(axis=0)
    sig = floor_r - P_prime
    sig_neg = (sig < 0).astype(float)  # indicator 1_{σ<0}

    # Compute matrices iteratively
    print(f"  Computing [{a},{b}] 4W={fourW} N={N}...", flush=True)
    G = np.zeros((N, N))
    K_coef = np.zeros((N, N))

    for j in range(N):
        pj = fej(RGRID - CEN[j]) + fej(RGRID + CEN[j])
        for k in range(N):
            pk = fej(RGRID - CEN[k]) + fej(RGRID + CEN[k])
            G[j, k] = (pj * pk).sum() * dr
            K_coef[j, k] = (pj * sig_neg * pk).sum() * dr

    # Solve generalized eigenvalue problem: K_coef @ c = λ G @ c
    try:
        eigs_val, eigs_vec = eigh(K_coef, G)
        # Sort descending
        idx = np.argsort(eigs_val)[::-1]
        eigs_val_desc = eigs_val[idx]

        lam_max = eigs_val_desc[0]
        min_sig = sig.min()
        floor = floor_r.mean()

        # Positivity certificate margin: floor - λ_max * |min σ|
        gap = floor - lam_max * np.abs(min_sig)

        # Window length L and time-bandwidth product Δt
        L = b - a
        delt = (L * 2*W) / np.pi

        return {
            'a': a, 'b': b, 'L': L, 'fourW': fourW, 'W': W, 'N': N, 'dr': dr,
            'lam_max': lam_max,
            'min_sigma': min_sig,
            'floor': floor,
            'gap': gap,
            'eigs_top_10': eigs_val_desc[:min(10, len(eigs_val_desc))],
            'neg_mass_max': lam_max * np.abs(min_sig),
            'neg_mass_frac_of_floor': (lam_max * np.abs(min_sig)) / floor,
            'delt_time_bandwidth': delt,
            'LW_plunge_region': 'critical' if delt < 1.5 else 'robust',
            'K_rank_numeric': np.sum(eigs_val_desc > 1e-9),
            'trace_K': K_coef.trace(),
            'norm_K': np.linalg.norm(K_coef),
        }
    except Exception as e:
        print(f"  ERROR at [{a},{b}] 4W={fourW} N={N}: {e}", flush=True)
        return None

def certification_check(result):
    """Check if the concentration margin is positive (certification of positivity)."""
    if result is None:
        return None
    gap = result['gap']
    certified = gap >= -1e-10  # allow tiny numerical error
    return {
        'window': f"[{result['a']},{result['b']}]",
        'fourW': result['fourW'],
        'N': result['N'],
        'floor': result['floor'],
        'lam_max': result['lam_max'],
        'min_sigma': result['min_sigma'],
        'neg_mass_max': result['neg_mass_max'],
        'gap': gap,
        'certified': certified,
        'delt': result['delt_time_bandwidth'],
        'LW_region': result['LW_plunge_region'],
    }

print("="*120, flush=True)
print("§888 CONCENTRATION OPERATOR K = B_W · 1_{σ<0} · B_W (in Fejér coefficient space)", flush=True)
print("λ_max(K) = max fraction of |φ|² on {σ<0}; positivity certificate: gap = floor − λ_max·|min σ| ⪰ 0", flush=True)
print("="*120, flush=True)

# Corner windows [14,44] with various 4W values
results_corner = []
for fourW in [2.0, 3.0, 4.0, 5.0]:
    for N in [12, 15, 17, 20]:
        # Ncap = L · 2W/π; keep N below Ncap
        Ncap = (44 - 14) * (fourW/4 * 2) / np.pi
        if N > Ncap * 1.1:  # allow 10% overshoot tolerance
            continue
        res = run_concentration(14, 44, fourW, N, dr=0.003)
        if res:
            results_corner.append(res)
            cert = certification_check(res)
            print(f"[14,44] 4W={fourW:4.1f} N={N:2d} (Ncap={Ncap:5.1f}): λ_max={res['lam_max']:.4f}  " +
                  f"min σ={res['min_sigma']:+.4f}  floor={res['floor']:.4f}  gap={cert['gap']:+.3e}  " +
                  f"CERT={cert['certified']}  Δt={res['delt_time_bandwidth']:.2f}")

print("\n" + "="*120)
print("Robust region windows (height climbs, LW plunge should flatten)")
print("="*120)

results_robust = []
for (a, b, fourW) in [(60, 90, 4.0), (120, 150, 4.0), (300, 330, 4.0), (1000, 1030, 4.0)]:
    N = max(20, int((b-a) * fourW/4 * 2 / np.pi * 0.8))  # heuristic N
    Ncap = (b - a) * (fourW/4 * 2) / np.pi
    if N > Ncap:
        N = int(Ncap * 0.9)
    if N < 8:
        N = 8

    res = run_concentration(a, b, fourW, N, dr=0.003)
    if res:
        results_robust.append(res)
        cert = certification_check(res)
        print(f"[{a:4d},{b:4d}] 4W={fourW:4.1f} N={N:2d} (Ncap={Ncap:5.1f}): λ_max={res['lam_max']:.4f}  " +
              f"min σ={res['min_sigma']:+.4f}  floor={res['floor']:.4f}  gap={cert['gap']:+.3e}  " +
              f"CERT={cert['certified']}  Δt={res['delt_time_bandwidth']:.2f}")

print("\n" + "="*120)
print("SUMMARY: CONCENTRATION OPERATOR ANALYSIS")
print("="*120)

# Compile all results
all_results = results_corner + results_robust
all_certs = [certification_check(r) for r in all_results]

print("\nCORNER [14,44] TRAJECTORY (4W varied):")
corner_4w_sorted = sorted([c for c in all_certs if c['window'] == '[14,44]'], key=lambda x: x['fourW'])
for c in corner_4w_sorted:
    print(f"  4W={c['fourW']:4.1f}: λ_max={c['lam_max']:.4f}  gap={c['gap']:+.3e}  Δt={c['delt']:.2f}  {c['certified']}")

print("\nHEIGHT SWEEP (4W=4.0 constant):")
height_4w4 = sorted([c for c in all_certs if c['fourW'] == 4.0], key=lambda x: int(c['window'].split(',')[0][1:]))
for c in height_4w4:
    print(f"  {c['window']:>10s}: λ_max={c['lam_max']:.4f}  gap={c['gap']:+.3e}  Δt={c['delt']:.2f}  {c['certified']}")

print("\n" + "="*120)
print("INTERPRETATION:")
print("="*120)

# Check for Landau-Widom plunge signature
lw_scores = []
for c in all_certs:
    if c['LW_region'] == 'critical':
        lw_scores.append(('critical', c['lam_max'], c['delt']))
    else:
        lw_scores.append(('robust', c['lam_max'], c['delt']))

critical_lams = [x[1] for x in lw_scores if x[0] == 'critical']
robust_lams = [x[1] for x in lw_scores if x[0] == 'robust']

if critical_lams and robust_lams:
    crit_mean = np.mean(critical_lams)
    rob_mean = np.mean(robust_lams)
    print(f"Mean λ_max at corner (Δt<1.5):   {crit_mean:.4f}")
    print(f"Mean λ_max in robust (Δt≥1.5):   {rob_mean:.4f}")
    print(f"Plunge ratio (critical/robust):  {crit_mean/rob_mean:.3f}")
    if abs(crit_mean - 0.5) < 0.15:
        print("✓ λ_max ≈ 0.5 at critical time-bandwidth product (LW plunge signature PRESENT)")
    else:
        print(f"  λ_max at critical: {crit_mean:.4f} (expected ~0.5 for LW plunge)")

print("\nCERTIFICATION STATUS:")
certified_count = sum(1 for c in all_certs if c['certified'])
total_count = len(all_certs)
print(f"  Positivity certificate (gap≥0): {certified_count}/{total_count} cases")

if certified_count == total_count:
    print("✓ FULL CERTIFICATION: all windows pass λ_max·|min σ| < floor")
    print("  → Positivity certificate margin ≥ 0 everywhere tested")
    print("  → Sufficient to prove RH corner (if Fejér completeness holds beyond 4W)")
elif certified_count > total_count * 0.8:
    print("⚠ PARTIAL CERTIFICATION: most windows pass, corner remains tight")
    print("  → Margin is positive but razor-thin at [14,44] corner")
else:
    print("✗ FAILED: negative gaps detected")
    print("  → This bound (Fejér K alone) is too weak to certify positivity")
    print("  → Need refined concentration bound or coupling-aware projection")

# Save detailed results
results_file = '/tmp/concentration_operator_results.json'
output_dict = {
    'results': [
        {
            'window': f"[{r['a']},{r['b']}]",
            'fourW': r['fourW'],
            'N': r['N'],
            'L': r['L'],
            'W': r['W'],
            'lam_max': float(r['lam_max']),
            'min_sigma': float(r['min_sigma']),
            'floor': float(r['floor']),
            'gap': float(r['gap']),
            'neg_mass_max': float(r['neg_mass_max']),
            'neg_mass_frac_of_floor': float(r['neg_mass_frac_of_floor']),
            'delt_time_bandwidth': float(r['delt_time_bandwidth']),
            'eigs_top_5': [float(x) for x in r['eigs_top_10'][:5]],
        }
        for r in all_results
    ],
    'certification_summary': {
        'total_cases': total_count,
        'certified_cases': certified_count,
        'all_positive': certified_count == total_count,
    }
}

with open(results_file, 'w') as f:
    json.dump(output_dict, f, indent=2)
print(f"\nDetailed results saved to {results_file}")
print("="*120)
