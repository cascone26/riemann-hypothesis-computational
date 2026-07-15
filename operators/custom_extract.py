#!/usr/bin/env python3
"""
custom_extract.py — Newton-based zero extraction from CvS eigenvector.

Uses our own F(z) implementation to avoid the extract_zeros sinc-node
proximity bug seen at k=9 (c=13,c=23) and predicted at k=7 (c=100).

Loads eigenvector from cvs_highprec_v2_c{c}.json (which saves eigvec).
Falls back to rebuilding the Galerkin matrix if eigenvec not saved.

Usage:
  python3 custom_extract.py <c> [T=800] [dps=150] [n_zeros=10]

Algorithm:
  1. Compute F(z) = sqrt(L) * sum_n c_n * exp(i*beta_n) * sinc_unnorm(beta_n)
     where beta_n = (2*pi*n - z*L) / 2
  2. Bracket minima of |F|^2 on a dense grid
  3. Refine with Newton on complex F(z) = 0
  4. Verify: report -log10(|F_detected|/|F_peak|) as extraction quality
"""
import sys, json, os, time
import mpmath as mp

# Parse args
args = sys.argv[1:]
c       = int(args[0])   if len(args) > 0 else 100
T       = int(args[1])   if len(args) > 1 else 800
dps     = int(args[2])   if len(args) > 2 else 150
n_zeros = int(args[3])   if len(args) > 3 else 10

mp.mp.dps = dps + 20
L = mp.log(c)
L_f = float(L)
TWO_PI = 2.0 * mp.pi

KNOWN_ZEROS_STR = [
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
KNOWN_ZEROS = [mp.mpf(s) for s in KNOWN_ZEROS_STR]

INFILE_V2 = os.path.expanduser(f"~/riemann/results/cvs_highprec_v2_c{c}.json")
OUTFILE   = os.path.expanduser(f"~/riemann/results/custom_extract_c{c}.json")


def load_eigvec():
    """Load eigenvector from v2 JSON (if saved)."""
    if not os.path.exists(INFILE_V2):
        return None
    print(f"Loading eigenvector from {INFILE_V2}...", flush=True)
    with open(INFILE_V2) as f:
        data = json.load(f)
    if 'eigvec' not in data:
        return None
    eigvec_strs = data['eigvec']
    # eigvec is list of (real, imag) string pairs
    N = (len(eigvec_strs) - 1) // 2
    idx = list(range(-N, N+1))
    eigvec = [mp.mpc(mp.mpf(re), mp.mpf(im)) for (re, im) in eigvec_strs]
    return idx, eigvec


def build_eigvec_from_scratch():
    """Rebuild Galerkin matrix and compute eigenvector."""
    try:
        from connes_cvs import build_galerkin_matrix, compute_ground_state
    except ImportError:
        print("ERROR: connes-cvs not installed and no saved eigvec.", file=sys.stderr)
        sys.exit(1)
    print(f"Rebuilding Galerkin matrix c={c} N=100 T={T} dps={dps}...", flush=True)
    t0 = time.time()
    Q = build_galerkin_matrix(c=c, N=100, T=T, dps=dps)
    print(f"  Build: {time.time()-t0:.1f}s", flush=True)
    lam_min, eigvec = compute_ground_state(Q)
    print(f"  lam_min = {mp.nstr(lam_min, 6)}", flush=True)
    N = 100
    idx = list(range(-N, N+1))
    return idx, list(eigvec)


def F_eval(z, idx, eigvec):
    """Compute F(z) using mpmath. F is approximately zero at Riemann zeros."""
    z = mp.mpc(z)
    s = mp.mpc(0)
    L_mp = mp.log(c)
    for n, cn in zip(idx, eigvec):
        beta = (2 * mp.pi * n - z * L_mp) / 2
        if mp.fabs(beta) < 1e-15:
            f = mp.mpc(1, 0)  # sinc(0) = 1, exp(0) = 1
        else:
            f = mp.exp(1j * beta) * mp.sin(beta) / beta
        s += cn * f
    return mp.sqrt(L_mp) * s


def dF_dz_eval(z, idx, eigvec):
    """Compute dF/dz."""
    z = mp.mpc(z)
    L_mp = mp.log(c)
    s = mp.mpc(0)
    for n, cn in zip(idx, eigvec):
        beta = (2 * mp.pi * n - z * L_mp) / 2
        dbeta_dz = -L_mp / 2
        if mp.fabs(beta) < 1e-15:
            # d/dbeta[e^{i*beta} * sin(beta)/beta] at beta→0 = i (from Taylor series)
            df_dbeta = mp.mpc(0, 1)
        else:
            # d/dbeta[e^{ibeta} * sin(beta)/beta]
            # = i*e^{ibeta}*sin(beta)/beta + e^{ibeta}*(cos(beta)/beta - sin(beta)/beta^2)
            eib = mp.exp(1j * beta)
            df_dbeta = 1j * eib * mp.sin(beta) / beta + eib * (mp.cos(beta) / beta - mp.sin(beta) / beta**2)
        s += cn * df_dbeta * dbeta_dz
    return mp.sqrt(L_mp) * s


def grid_search_bracket(idx, eigvec, z_lo, z_hi, n_pts=3000):
    """Dense grid search for |F| minima as initial zero candidates."""
    import numpy as np
    z_grid = np.linspace(float(z_lo), float(z_hi), n_pts)
    # Fast double-precision evaluation for bracketing
    L_f = float(mp.log(c))
    n_arr = np.array([float(n) for n in idx])
    c_arr_re = np.array([float(mp.re(cn)) for cn in eigvec])
    c_arr_im = np.array([float(mp.im(cn)) for cn in eigvec])
    c_arr = c_arr_re + 1j * c_arr_im

    absF = np.zeros(n_pts)
    for i, z in enumerate(z_grid):
        beta = (2 * np.pi * n_arr - z * L_f) / 2
        with_eps = np.where(np.abs(beta) < 1e-14, 1e-14, beta)
        f = np.exp(1j * beta) * np.sin(with_eps) / with_eps
        F = np.sqrt(L_f) * np.dot(c_arr, f)
        absF[i] = abs(F)

    # Find local minima
    mins = []
    for i in range(1, n_pts-1):
        if absF[i] < absF[i-1] and absF[i] < absF[i+1]:
            mins.append((absF[i], z_grid[i]))
    mins.sort()
    return mins


def newton_refine(z0, idx, eigvec, max_iter=60, tol=None):
    """Newton's method on F(z)=0 starting from z0."""
    if tol is None:
        tol = mp.mpf(10) ** (-(dps - 5))
    z = mp.mpc(z0)
    for it in range(max_iter):
        Fz = F_eval(z, idx, eigvec)
        dFz = dF_dz_eval(z, idx, eigvec)
        if mp.fabs(dFz) < 1e-300:
            break
        step = Fz / dFz
        z = z - step
        if mp.fabs(step) < tol:
            break
    return mp.re(z)  # zeros on real axis, take real part


def extract_zeros_custom(idx, eigvec, n_zeros, z_lo=12.0, z_hi=55.0):
    """Extract first n_zeros Riemann zeros via grid search + Newton."""
    print(f"\nGrid search [{z_lo}, {z_hi}]...", flush=True)
    mins = grid_search_bracket(idx, eigvec, z_lo, z_hi, n_pts=5000)

    # Filter: keep candidates within ±3 of known zeros
    results = []
    used = set()
    for absval, z_cand in mins:
        # Find nearest known zero
        for k_idx, gz in enumerate(KNOWN_ZEROS[:n_zeros]):
            gz_f = float(gz)
            if abs(z_cand - gz_f) < 2.0 and k_idx not in used:
                used.add(k_idx)
                print(f"  k={k_idx+1}: initial z={z_cand:.6f} (|F|={absval:.3e}), refining...", flush=True)
                z_refined = newton_refine(z_cand, idx, eigvec)
                Fz_refined = F_eval(z_refined, idx, eigvec)
                results.append((k_idx+1, z_refined, Fz_refined))
                break

    # Sort by k
    results.sort(key=lambda x: x[0])
    return results


def main():
    print(f"\ncustom_extract: c={c}, T={T}, dps={dps}, n_zeros={n_zeros}")
    print(f"L = log({c}) = {mp.nstr(mp.log(c), 10)}")

    # Load or build eigenvector
    loaded = load_eigvec()
    if loaded is None:
        print("No saved eigvec found, rebuilding from scratch...")
        idx, eigvec = build_eigvec_from_scratch()
    else:
        idx, eigvec = loaded

    print(f"Eigenvector: {len(eigvec)} components (idx {idx[0]} to {idx[-1]})", flush=True)

    # Compute max |F| for normalization
    import numpy as np
    L_f = float(mp.log(c))
    n_arr = np.array([float(n) for n in idx])
    c_arr = np.array([float(mp.re(cn)) + 1j * float(mp.im(cn)) for cn in eigvec])
    z_norm = np.linspace(12, 55, 1000)
    absF_max = 0.0
    for z in z_norm:
        beta = (2 * np.pi * n_arr - z * L_f) / 2
        with_eps = np.where(np.abs(beta) < 1e-14, 1e-14, beta)
        F_dp = np.sqrt(L_f) * np.dot(c_arr, np.exp(1j * beta) * np.sin(with_eps) / with_eps)
        absF_max = max(absF_max, abs(F_dp))

    # Extract zeros
    t0 = time.time()
    zero_list = extract_zeros_custom(idx, eigvec, n_zeros)
    dt = time.time() - t0

    # Report
    print(f"\n{'k':>4}  {'gamma_detected':>50}  {'error':>20}  {'digits':>8}  {'|F|/peak':>12}")
    zero_results = []
    for k, z_det, Fz_det in zero_list:
        if k <= len(KNOWN_ZEROS):
            gz = KNOWN_ZEROS[k-1]
            err = z_det - gz
            digits = -float(mp.log10(mp.fabs(err / gz))) if mp.fabs(err) > 0 else dps
            # |F|/peak
            Fz_abs_f = float(mp.fabs(Fz_det))
            depth = Fz_abs_f / absF_max
            print(f"{k:>4}  {mp.nstr(z_det, min(50, dps)):<50}  {mp.nstr(err, 4):>20}  {digits:>8.1f}  {depth:>12.3e}")
            zero_results.append({
                'k': k,
                'gamma_detected': mp.nstr(z_det, dps),
                'gamma_true': KNOWN_ZEROS_STR[k-1],
                'error': mp.nstr(err, 4),
                'digits': digits,
                'F_depth': depth,
            })
        else:
            print(f"{k:>4}  {mp.nstr(z_det, 20):<50}")

    out = {
        'c': c, 'T': T, 'dps': dps, 'n_zeros': n_zeros,
        'method': 'custom_newton',
        'dt_extract': round(dt, 1),
        'zeros': zero_results,
    }
    with open(OUTFILE, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUTFILE} (extract time: {dt:.1f}s)")


if __name__ == "__main__":
    main()
