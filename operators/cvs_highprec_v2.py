#!/usr/bin/env python3
"""
cvs_highprec_v2.py — High-precision CvS extraction, saves eigenvector for custom extraction.

Same as cvs_highprec_olympus.py but also saves the eigenvector to JSON
so custom_extract.py can re-run zero extraction without rebuilding the matrix.

Usage: python3 cvs_highprec_v2.py [c] [N] [T] [dps] [n_zeros]
Defaults: c=100, N=100, T=800, dps=150, n_zeros=10

Saves:
  ~/riemann/results/cvs_highprec_v2_c{c}.json  (includes eigvec)
"""
import sys, json, os, time
import mpmath as mp

try:
    from connes_cvs import build_galerkin_matrix, compute_ground_state, extract_zeros
except ImportError:
    print("ERROR: connes-cvs not installed.")
    sys.exit(1)

args = sys.argv[1:]
c       = int(args[0]) if len(args) > 0 else 100
N       = int(args[1]) if len(args) > 1 else 100
T       = int(args[2]) if len(args) > 2 else 800
dps     = int(args[3]) if len(args) > 3 else 150
n_zeros = int(args[4]) if len(args) > 4 else 10

KNOWN_ZEROS = [
    "14.134725141734693790457251983562470270784257115699243175685567460149963429809256764949010393171561012779202971548797484545370608645013143",
    "21.022039638771554992628479593896902777334340524902781754629520403501896418518816088404191854501501765040284052568239396498454078839951492",
    "25.010857580145688763213790992562821818659549672557996672496542006745092098441644890444124638720514698760252982851416978174939717826551490",
    "30.424876125859513210311897530584091320181560023715440180962146036993329389955868978867882906484166760283875069015936069578844648882476729",
    "32.935061587739189690662368964074903488812715603517039009280003440784815608630312756069569059028778703483769476143097023048834707723219543",
    "37.586178158825671257217763480705332821405597350830793218333001113622369333985930842675015390924628988399942285372892798810696185639285508",
    "40.918719012147495187398126914633254395726165962777279536161303667253767216537228457284000560997774727938018506268988285745949843958218726",
    "43.327073280914999519496122165406805782645668371836871446099472580924007635987979516699888940007520261862068668768773413700990997218734283",
    "48.005150881167159727942472749427516588896474580064948426499037543597855756568560157905047949046047034831048451393348572723893948046703553",
    "49.773832477672302181916784678563724057723178299676662796994823491174419282791369879682893278028684428838823388001782063044950827427282929",
]

OUTFILE = os.path.expanduser(f"~/riemann/results/cvs_highprec_v2_c{c}.json")
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

print(f"cvs_highprec_v2: c={c}, N={N}, T={T}, dps={dps}, n_zeros={n_zeros}")
print(f"Output: {OUTFILE}", flush=True)

mp.mp.dps = dps + 20
L = mp.log(c)
print(f"L = {mp.nstr(L, 8)}", flush=True)

t0 = time.time()
print("Building Galerkin matrix...", flush=True)
Q = build_galerkin_matrix(c=c, N=N, T=T, dps=dps)
t1 = time.time()
print(f"  Done in {t1-t0:.1f}s. Shape: {Q.rows}x{Q.cols}", flush=True)

print("Computing ground state...", flush=True)
lam_min, eigvec = compute_ground_state(Q)
t2 = time.time()
lam_digits = -float(mp.log10(mp.fabs(lam_min)))
print(f"  lam_min = {mp.nstr(lam_min, 6)}, digits = {lam_digits:.1f} ({t2-t1:.1f}s)", flush=True)

print(f"\nExtracting {n_zeros} zeros...", flush=True)
t3 = time.time()
extracted = extract_zeros(eigvec, L=L, n_zeros=n_zeros, dps=dps)
t4 = time.time()
print(f"  Done in {t4-t3:.1f}s", flush=True)

print(f"\n{'k':>4}  {'gamma':>50}  {'error':>20}  {'digits':>8}")
zero_results = []
for ze in extracted:
    k  = ze['k']
    gd = ze.get('gamma_detected')
    if k <= len(KNOWN_ZEROS) and gd is not None:
        gk_true = mp.mpf(KNOWN_ZEROS[k-1])
        err = gd - gk_true
        digits = -float(mp.log10(mp.fabs(err / gk_true))) if mp.fabs(err) > 0 else dps
        print(f"{k:>4}  {mp.nstr(gd, min(50, dps)):<50}  {mp.nstr(err, 4):>20}  {digits:>8.1f}")
        zero_results.append({
            'k': k,
            'gamma_detected': mp.nstr(gd, dps),
            'gamma_true': KNOWN_ZEROS[k-1],
            'error': mp.nstr(err, 4),
            'digits': digits,
        })
    else:
        print(f"{k:>4}  {mp.nstr(gd, 15) if gd else 'None':>50}")

# Save eigenvector as (re, im) string pairs for each component
print("\nSerializing eigenvector...", flush=True)
eigvec_list = [[mp.nstr(mp.re(x), dps), mp.nstr(mp.im(x), dps)] for x in eigvec]

out = {
    'c': c, 'N': N, 'T': T, 'dps': dps, 'n_zeros': n_zeros,
    'L': mp.nstr(L, 10),
    'lam_min': mp.nstr(lam_min, 8),
    'lam_digits': lam_digits,
    'dt_build': round(t1 - t0, 1),
    'dt_eig': round(t2 - t1, 1),
    'dt_zeros': round(t4 - t3, 1),
    'zeros': zero_results,
    'eigvec': eigvec_list,
}

with open(OUTFILE, 'w') as f:
    json.dump(out, f, indent=2)

# Also save without eigvec for quick inspection
out_no_eigvec = {k: v for k, v in out.items() if k != 'eigvec'}
with open(OUTFILE.replace('.json', '_summary.json'), 'w') as f:
    json.dump(out_no_eigvec, f, indent=2)

print(f"\nSaved: {OUTFILE}")
print(f"Saved: {OUTFILE.replace('.json','_summary.json')} (no eigvec)")
