#!/usr/bin/env python3
"""ANGLE A REFINED: Measure VANISHING more carefully. Instead of discrete grid evaluation,
measure the decay of |φ(r)| near the zeros, vs far from zeros. Also check whether the eigenvector
coefficients themselves show concentration/sparsity patterns that relate to zeros."""
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

def analyze_corner_v2(a, b, fourW, frac=0.65, dr=0.004):
    """Build M_σ, whiten, eigendecompose; analyze near-null modes with multiple vanishing metrics."""
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

    # Near-null cluster (positive, small)
    nn_mask = (ev > 1e-8) & (ev < 1e-4)
    nn_indices = np.where(nn_mask)[0]
    nn_eigenvalues = ev[nn_indices]
    nn_eigenvectors = V[:, nn_indices]

    zeta_zeros = zeros_in(a, b)
    num_zeros = len(zeta_zeros)

    # Unwhiten
    Gsqrt = gv @ np.diag(np.sqrt(gw)) @ gv.T

    # Reconstruct band-limited functions
    phi_vals_full = {}  # eigvec index -> full φ on RGRID
    phi_at_zeros = {}
    phi_at_far = {}
    vanishing_data = []

    for i, (eig_val, v_w) in enumerate(zip(nn_eigenvalues, nn_eigenvectors.T)):
        v_orig = Gsqrt @ v_w
        phi_full = (v_orig[:, None] * P).sum(axis=0)  # shape (len(RGRID),)
        phi_vals_full[i] = phi_full

        # Evaluation near zeros vs far
        # For each zero, find nearest grid points within ±W (bandwidth scale)
        near_window = W * 2  # ±2 * bandwidth
        phi_at_zero_list = []
        phi_near_list = []

        for z in zeta_zeros:
            # Closest grid point
            idx_closest = np.argmin(np.abs(RGRID - z))
            phi_at_zero_list.append(phi_full[idx_closest])

            # Points within ±2*bandwidth
            mask_near = np.abs(RGRID - z) < near_window
            if mask_near.sum() > 0:
                phi_near_list.extend(phi_full[mask_near])

        # Decay metric: RMS near zero vs RMS far from any zero
        if len(phi_at_zero_list) > 0:
            rms_at_zeros = np.sqrt(np.mean(np.array(phi_at_zero_list)**2))
        else:
            rms_at_zeros = 0

        # Far from zeros: points that are >4*bandwidth away from all zeros
        far_window = 4 * W
        mask_far = np.ones(len(RGRID), dtype=bool)
        for z in zeta_zeros:
            mask_far &= np.abs(RGRID - z) > far_window
        if mask_far.sum() > 0:
            rms_far = np.sqrt(np.mean(phi_full[mask_far]**2))
        else:
            rms_far = np.sqrt(np.mean(phi_full**2))

        decay_ratio = rms_at_zeros / rms_far if rms_far > 0 else 0

        # Eigenvector sparsity: how many basis functions (CEN) have large coefficients?
        v_abs = np.abs(v_orig)
        threshold = np.median(v_abs)  # median threshold
        num_large = (v_abs > threshold).sum()
        sparsity = num_large / N

        vanishing_data.append({
            'eig_val': eig_val,
            'rms_at_zeros': rms_at_zeros,
            'rms_far': rms_far,
            'decay_ratio': decay_ratio,
            'eigvec_sparsity': sparsity,
            'num_large_coeff': num_large
        })

    return {
        'corner': (a, b),
        'fourW': fourW,
        'N': N,
        'Ncap': Ncap,
        'num_zeros': num_zeros,
        'num_near_null': len(nn_eigenvalues),
        'near_null_eigenvalues': nn_eigenvalues,
        'vanishing_data': vanishing_data,
        'zeta_zeros': zeta_zeros,
        'CEN': CEN
    }

corners = [
    (14, 44, 4.0),
    (14, 34, 4.0),
    (14, 44, 5.0)
]

print("="*140)
print("ANGLE A REFINED: NEAR-NULL VANISHING — DECAY METRIC & SPARSITY")
print("Q: Do near-null eigenvectors decay at ζ-zeros and show sparse structure?")
print("="*140)

all_results = []
for a, b, fourW in corners:
    result = analyze_corner_v2(a, b, fourW)
    all_results.append(result)

    print(f"\n{'='*140}")
    print(f"CORNER: [{a}, {b}], 4W={fourW}")
    print(f"{'='*140}")
    print(f"Capacity N_cap = {result['Ncap']:.1f};  Actual N = {result['N']}")
    print(f"Number of ζ-zeros in window: {result['num_zeros']}")
    print(f"Number of near-null modes (λ ∈ [1e-8, 1e-4]): {result['num_near_null']}")

    if result['num_near_null'] > 0:
        print(f"\nDECAY ANALYSIS: RMS(φ at zeros) / RMS(φ far from zeros)")
        print(f"  {'idx':>3} {'λ':>12} {'RMS@zeros':>12} {'RMS_far':>12} {'ratio':>10} {'sparsity':>10} {'#large':>5}")
        for i, data in enumerate(result['vanishing_data']):
            print(f"  {i:>3} {data['eig_val']:>12.4e} {data['rms_at_zeros']:>12.4e} {data['rms_far']:>12.4e} {data['decay_ratio']:>10.4f} {data['eigvec_sparsity']:>10.4f} {data['num_large_coeff']:>5}")

        decay_ratios = np.array([d['decay_ratio'] for d in result['vanishing_data']])
        sparsities = np.array([d['eigvec_sparsity'] for d in result['vanishing_data']])

        print(f"\nSUMMARY STATS:")
        print(f"  Decay ratio (RMS@zeros/RMS_far):")
        print(f"    Mean: {decay_ratios.mean():.4f}")
        print(f"    Min:  {decay_ratios.min():.4f}")
        print(f"  Eigenvector sparsity (fraction with large coeffs):")
        print(f"    Mean: {sparsities.mean():.4f}")
        print(f"    Min:  {sparsities.min():.4f}")

        # Test hypothesis: if decay_ratio << 1, the modes vanish at zeros
        strong_decay = (decay_ratios < 0.1).sum()
        print(f"\n  Hypothesis Test: decay_ratio << 1 (e.g., < 0.1)?")
        print(f"    {strong_decay} / {len(decay_ratios)} modes with ratio < 0.1 → ", end="")
        if strong_decay / len(decay_ratios) > 0.5:
            print("YES (>50%) — modes decay at zeros")
        else:
            print("NO (<50%) — modes do NOT decay at zeros")
    else:
        print("  [No near-null modes]")

print(f"\n{'='*140}")
print("SUMMARY: Does #near-null ≈ #zeros AND decay_ratio << 1?")
print(f"{'='*140}")
print(f"{'corner':>12} {'4W':>4} {'#zeros':>7} {'#NN':>4} {'ratio':>8} {'mean_decay':>12}")
for result in all_results:
    a, b = result['corner']
    mean_decay = np.mean([d['decay_ratio'] for d in result['vanishing_data']]) if result['num_near_null'] > 0 else np.nan
    ratio_str = f"{result['num_near_null']}/{result['num_zeros']}"
    print(f"{f'[{a},{b}]':>12} {result['fourW']:>4.1f} {result['num_zeros']:>7} {result['num_near_null']:>4} {ratio_str:>8} {mean_decay:>12.4f}")

print(f"{'='*140}")
print("If decay_ratio ~1.0 (or >0.5), the hypothesis FAILS: near-null modes do NOT vanish at ζ-zeros.")
print("This would suggest the §855 arithmetic identity requires refinement or different interpretation.")
print(f"{'='*140}")
