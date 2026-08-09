#!/usr/bin/env python3
"""ANGLE A FINAL: Test whether near-null eigenvector COEFFICIENTS (not reconstructed functions)
show structure correlated with ζ-zero positions. If modes are "best vanishing" at zeros, perhaps
the basis function indices with large coefficients align with center positions BETWEEN zeros."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
import mpmath as mp

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam

NMAX = 4000; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

mp.mp.dps = 20
ZT = []
k = 1
while True:
    g = float(mp.im(mp.zetazero(k))); ZT.append(g)
    if g > 360: break
    k += 1
ZT = np.array(ZT)

def zeros_in(a, b):
    return ZT[(ZT >= a) & (ZT <= b)]

def analyze_corner_v3(a, b, fourW, frac=0.65, dr=0.004):
    """Analyze eigenvector coefficient structure vs ζ-zero positions."""
    W = fourW/4
    L = b - a
    Ncap = L * 2 * W / np.pi
    N = max(4, int(round(frac * Ncap)))

    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    CEN = np.linspace(a, b, N)

    def fej(x):
        y = W*x
        out = np.ones_like(y)
        nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2
        return out

    P = np.array([fej(RGRID-c) + fej(RGRID+c) for c in CEN])
    G = (P @ P.T) * dr
    arch = (P * OMEGA[None, :]) @ P.T * dr / (2*np.pi)

    mask = LOGN < fourW
    ln = LOGN[mask]
    lw = LAMW[mask]
    cosu = np.cos(np.outer(ln, RGRID))
    prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :] * P) @ cosu.T) * dr / (2*np.pi)
        prime[j, :] = 2.0 * (gp * lw[None, :]).sum(axis=1)
    prime = 0.5 * (prime + prime.T)

    def fc(c):
        y = W*(0.5j - c)
        y2 = W*(0.5j + c)
        return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2

    fp = np.array([fc(c) for c in CEN])
    pole = 2.0 * np.real(np.outer(fp, fp))

    gw, gv = eigh(G)
    B = gv @ np.diag(1/np.sqrt(gw)) @ gv.T
    Msig = arch - prime
    Msig_w = B @ Msig @ B
    ev, V = eigh(Msig_w)

    # Near-null cluster
    nn_mask = (ev > 1e-8) & (ev < 1e-4)
    nn_indices = np.where(nn_mask)[0]
    nn_eigenvalues = ev[nn_indices]
    nn_eigenvectors_w = V[:, nn_indices]  # in whitened coords

    # Unwhiten to original basis
    Gsqrt = gv @ np.diag(np.sqrt(gw)) @ gv.T
    nn_eigenvectors = np.array([Gsqrt @ v for v in nn_eigenvectors_w.T]).T  # shape (N, #nn)

    zeta_zeros = zeros_in(a, b)

    coeff_analysis = []
    for i, (eig_val, v_coeff) in enumerate(zip(nn_eigenvalues, nn_eigenvectors.T)):
        # Which basis functions (CEN indices) have the largest coefficients?
        abs_coeff = np.abs(v_coeff)
        sorted_idx = np.argsort(-abs_coeff)  # descending
        top_k = max(1, N // 2)  # top 50%
        top_indices = sorted_idx[:top_k]
        top_centers = CEN[top_indices]

        # How far are the top centers from the nearest ζ-zero?
        if len(zeta_zeros) > 0:
            distances_to_zeros = []
            for c in top_centers:
                min_dist = np.min(np.abs(c - zeta_zeros))
                distances_to_zeros.append(min_dist)
            distances_to_zeros = np.array(distances_to_zeros)
            mean_dist = distances_to_zeros.mean()
            min_dist_overall = distances_to_zeros.min()

            # Also check: are top centers BETWEEN zeros or AT zeros?
            # Define "between" as closer to midpoint of two adjacent zeros than to either zero
            between_count = 0
            for c in top_centers:
                # Find nearest two zeros
                idx_near = np.argsort(np.abs(c - zeta_zeros))[:2]
                if len(idx_near) >= 2:
                    z1, z2 = zeta_zeros[idx_near[0]], zeta_zeros[idx_near[1]]
                    midpoint = (z1 + z2) / 2
                    if np.abs(c - midpoint) < np.abs(c - z1):
                        between_count += 1

            coeff_analysis.append({
                'eig_val': eig_val,
                'num_large': top_k,
                'top_centers': top_centers,
                'mean_dist_to_zero': mean_dist,
                'min_dist_to_zero': min_dist_overall,
                'between_count': between_count,
                'between_ratio': between_count / top_k if top_k > 0 else 0
            })
        else:
            coeff_analysis.append({
                'eig_val': eig_val,
                'num_large': top_k,
                'top_centers': top_centers,
                'mean_dist_to_zero': np.nan,
                'min_dist_to_zero': np.nan,
                'between_count': 0,
                'between_ratio': np.nan
            })

    return {
        'corner': (a, b),
        'fourW': fourW,
        'N': N,
        'Ncap': Ncap,
        'num_zeros': len(zeta_zeros),
        'num_near_null': len(nn_eigenvalues),
        'near_null_eigenvalues': nn_eigenvalues,
        'coeff_analysis': coeff_analysis,
        'zeta_zeros': zeta_zeros,
        'CEN': CEN
    }

corners = [
    (14, 44, 4.0),
    (14, 34, 4.0),
    (14, 44, 5.0)
]

print("="*150)
print("ANGLE A FINAL: EIGENVECTOR COEFFICIENT STRUCTURE vs ζ-ZEROS")
print("Q: Do basis functions with large coefficients BETWEEN/AT zeros, suggesting vanishing structure?")
print("="*150)

all_results = []
for a, b, fourW in corners:
    result = analyze_corner_v3(a, b, fourW)
    all_results.append(result)

    print(f"\n{'='*150}")
    print(f"CORNER: [{a}, {b}], 4W={fourW}")
    print(f"{'='*150}")
    print(f"Capacity N_cap = {result['Ncap']:.1f};  Actual N = {result['N']}")
    print(f"Number of ζ-zeros in window: {result['num_zeros']}")
    print(f"Zeros at: {np.array2string(result['zeta_zeros'], precision=2)}")
    print(f"Basis centers (CEN): {np.array2string(result['CEN'], precision=2)}")
    print(f"Number of near-null modes: {result['num_near_null']}")

    if result['num_near_null'] > 0:
        print(f"\nCOEFFICIENT STRUCTURE (top 50% by magnitude):")
        print(f"  {'idx':>3} {'λ':>12} {'top_centers':>50} {'mean_dist_to_0':>15} {'between':>8}")
        for i, analysis in enumerate(result['coeff_analysis']):
            tc_str = np.array2string(analysis['top_centers'], precision=2, separator=', ')
            md_str = f"{analysis['mean_dist_to_zero']:.3f}" if not np.isnan(analysis['mean_dist_to_zero']) else "—"
            between_str = f"{analysis['between_ratio']:.2f}" if not np.isnan(analysis['between_ratio']) else "—"
            print(f"  {i:>3} {analysis['eig_val']:>12.4e} {tc_str:>50} {md_str:>15} {between_str:>8}")

        # Summary
        mean_dists = np.array([a['mean_dist_to_zero'] for a in result['coeff_analysis'] if not np.isnan(a['mean_dist_to_zero'])])
        between_ratios = np.array([a['between_ratio'] for a in result['coeff_analysis'] if not np.isnan(a['between_ratio'])])

        print(f"\nSTRUCTURE SUMMARY:")
        if len(mean_dists) > 0:
            print(f"  Mean distance of top centers to nearest zero: {mean_dists.mean():.3f}")
            print(f"  Typical zero spacing in window: {(result['zeta_zeros'].max() - result['zeta_zeros'].min()) / max(1, result['num_zeros']-1):.3f}")
        if len(between_ratios) > 0:
            print(f"  Fraction of top centers BETWEEN adjacent zeros: {between_ratios.mean():.3f}")
            print(f"    → Strong between-zero concentration? ", end="")
            if between_ratios.mean() > 0.6:
                print("YES — suggests modes optimize to vanish AT zeros")
            else:
                print("NO — modes not concentrated between zeros")

print(f"\n{'='*150}")
print("CONCLUSION: Do eigenvector coefficients show structure (between zeros, far from zeros)?")
print(f"{'='*150}")
print(f"{'corner':>12} {'4W':>4} {'#zeros':>7} {'#NN':>4} {'mean_dist':>12} {'between%':>10}")
for result in all_results:
    a, b = result['corner']
    if result['num_near_null'] > 0:
        mean_dists = np.array([a['mean_dist_to_zero'] for a in result['coeff_analysis'] if not np.isnan(a['mean_dist_to_zero'])])
        between_ratios = np.array([a['between_ratio'] for a in result['coeff_analysis'] if not np.isnan(a['between_ratio'])])
        md = mean_dists.mean() if len(mean_dists) > 0 else np.nan
        br = between_ratios.mean() if len(between_ratios) > 0 else np.nan
        print(f"{f'[{a},{b}]':>12} {result['fourW']:>4.1f} {result['num_zeros']:>7} {result['num_near_null']:>4} {md:>12.3f} {br*100:>9.1f}%")

print(f"{'='*150}")
print("If between% >> 50% or mean_dist ~0, the modes ARE structured around zeros (§855 ✓).")
print("If mean_dist ~ window_size/2 or between% ~0%, the modes are NOT structured (§855 ✗).")
print(f"{'='*150}")
