#!/usr/bin/env python3
"""ANGLE 2: Is M_σ's near-null cluster the Landau–Widom prolate plunge?
The near-null cluster of M_σ (eigenvalues ~1e-7) carries the RH razor.
Conjecture: this is the prolate-spheroidal / time-frequency concentration plunge.

TEST: (1) Build concentration operator C = P_band ∘ (multiply by 1_window) ∘ P_band.
     Its eigenvalues give the time-frequency concentration: cluster ~1 (well-concentrated),
     sharp plunge (transition), cluster ~0 (poorly concentrated).
     (2) Compare the plunge region with M_σ's near-null cluster:
         - dimension match?
         - does M_σ's negative eigvec overlap with the plunge eigenspace?
         - plunge width vs Landau–Widom estimate ln(TB)/π²?
     (3) Report prolate spectrum, plunge width, overlaps — HANDLE or negative result (REPORT).
"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
np.set_printoptions(precision=5, suppress=True, linewidth=140)

def vonmangoldt(N):
    """von Mangoldt function: Λ(n) = log p if n=p^k, else 0."""
    lam = np.zeros(N+1)
    isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True
            lp = np.log(p)
            pk = p
            while pk <= N:
                lam[pk] = lp
                pk *= p
    return lam

def run(a, b, fourW, frac=0.65, dr=0.004):
    """
    Build Fejér band-limited basis, compute M_σ (arch-prime),
    and construct the concentration operator C.

    Returns: (1) prolate spectrum (concentration eigenvalues),
             (2) plunge region characterization,
             (3) comparison with M_σ's near-null cluster.
    """
    W = fourW / 4
    L = b - a
    Ncap = L * 2 * W / np.pi
    N = max(4, int(round(frac * Ncap)))

    # Compute grid for Fejér band-limited functions
    RGRID = np.arange(-(b+60), (b+60), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    CEN = np.linspace(a, b, N)

    # Fejér kernel: sinc^2(W*x)
    def fej(x):
        y = W * x
        out = np.ones_like(y)
        nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz]) / y[nz]) ** 2
        return out

    # P[j,k] = fej(RGRID[k] - CEN[j]) + fej(RGRID[k] + CEN[j])
    # (includes both ±c to get symmetry)
    P = np.array([fej(RGRID - c) + fej(RGRID + c) for c in CEN])

    # Gram matrix G and archimedean form
    G = (P @ P.T) * dr
    arch = (P * OMEGA[None, :]) @ P.T * dr / (2 * np.pi)

    # Prime form: Σ_{log n < 4W} Λ(n)/√n cos(r log n)
    NMAX = 4000
    LAM = vonmangoldt(NMAX)
    NS = np.arange(2, NMAX+1)
    LOGN = np.log(NS)
    LAMW = LAM[2:] / np.sqrt(NS)

    mask = LOGN < fourW
    ln = LOGN[mask]
    lw = LAMW[mask]
    cosu = np.cos(np.outer(ln, RGRID))
    prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :] * P) @ cosu.T) * dr / (2 * np.pi)
        prime[j, :] = 2.0 * (gp * lw[None, :]).sum(axis=1)
    prime = 0.5 * (prime + prime.T)

    # M_σ = arch - prime
    Msig = arch - prime

    # === CONCENTRATION OPERATOR ===
    # Classical prolate spheroidal problem:
    # - C = "time-limited" operator: multiply by window indicator 1_{[a,b]}
    # - G = Gram matrix (normalization)
    # - Solve generalized eigenvalue problem: C v = λ G v
    # This gives concentration eigenvalues ∈ [0,1], where λ ≈ 1 means well-concentrated
    #   (function is mostly inside the window) and λ ≈ 0 means poorly concentrated.

    window_mask = (RGRID >= a) & (RGRID <= b)
    C = (P * window_mask[None, :]) @ P.T * dr
    C = 0.5 * (C + C.T)  # symmetrize

    # Solve generalized eigenvalue problem: C v = λ G v
    # Use eigh with the b parameter
    from scipy.linalg import eigh as scipy_eigh
    concen_evals_gen, concen_evecs_gen = scipy_eigh(C, G)
    concen_evals_gen = np.clip(concen_evals_gen, 0, 1)  # values should be in [0,1]

    # Gram whitening for fair comparison with M_σ
    gw, gv = eigh(G)
    Ginv_sqrt = gv @ np.diag(1.0 / np.sqrt(gw)) @ gv.T

    # === EXTRACT PLUNGE REGION ===
    # Sort eigenvalues and find the plunge (region where they drop sharply)
    sorted_idx = np.argsort(-concen_evals_gen)  # descending
    sorted_evals = concen_evals_gen[sorted_idx]
    sorted_evecs = concen_evecs_gen[:, sorted_idx]

    # Plunge region: find the largest spectral gap (steepest drop)
    # This is more robust than arbitrary thresholds
    if len(sorted_evals) > 1:
        gaps = np.diff(sorted_evals[::-1])  # ascending order
        max_gap_idx = np.argmax(gaps)
        # Plunge starts after the largest gap
        plunge_start = max_gap_idx + 1
        plunge_end = len(sorted_evals)
        plunge_indices = np.arange(plunge_start, plunge_end)
        plunge_mask = np.zeros(len(sorted_evals), dtype=bool)
        plunge_mask[plunge_indices] = True
        plunge_count = len(plunge_indices)
        plunge_max_gap = float(gaps[max_gap_idx])
    else:
        plunge_indices = np.array([])
        plunge_mask = np.zeros(len(sorted_evals), dtype=bool)
        plunge_count = 0
        plunge_max_gap = 0.0

    # Landau–Widom estimate for plunge width
    # For a window [a,b] with bandwidth W, the number of "prolate plunge" eigenvalues
    # is approximately ln(2TB) / π² or similar. Here T = L, B = 2W (since we cover ±W).
    T = L
    B = 2 * W
    lw_estimate = np.log(2 * T * B) / (np.pi ** 2)

    # === COMPARE WITH M_σ'S NEAR-NULL CLUSTER ===
    # M_σ's near-null cluster has eigenvalues ~1e-7
    # This is in the full basis, not G-normalized. We need to compare fairly.

    # Whiten M_σ to G-metric
    Msig_w = gv.T @ Msig @ gv / np.sqrt(np.outer(gw, gw))
    # Or more carefully: B @ M_σ @ B where B = G^{-1/2}
    B = Ginv_sqrt
    Msig_normalized = B @ Msig @ B

    ev_Msig, evec_Msig = eigh(Msig_normalized)

    # Near-null cluster of M_σ: eigenvalues close to 0
    # Extract the negative eigenvalue
    neg_mask_Msig = ev_Msig < -1e-10
    nneg = neg_mask_Msig.sum()

    if nneg > 0:
        mu = -float(ev_Msig[0])  # magnitude of negative eigenvalue
        w = evec_Msig[:, 0]  # eigenvector of negative eigenvalue
    else:
        mu = 0.0
        w = np.zeros(N)

    # Near-null cluster: find by spectral gap (similar to plunge detection)
    # The near-null cluster is the set of eigenvalues where there's a gap separating them
    # from the bulk of the spectrum
    sorted_Msig_idx = np.argsort(np.abs(ev_Msig))  # ascending by magnitude
    sorted_Msig_evals = ev_Msig[sorted_Msig_idx]
    sorted_Msig_evecs = evec_Msig[:, sorted_Msig_idx]

    if len(sorted_Msig_evals) > 1:
        abs_gaps = np.diff(np.abs(sorted_Msig_evals))
        max_gap_idx = np.argmax(abs_gaps)
        # Near-null cluster is before the largest gap
        near_null_end = max_gap_idx + 1
        near_null_indices = np.arange(0, near_null_end)
        near_null_count = near_null_end
    else:
        near_null_indices = np.array([0])
        near_null_count = 1

    near_null_evals = sorted_Msig_evals[:near_null_count]
    near_null_evecs = sorted_Msig_evecs[:, :near_null_count]

    # === OVERLAP TEST ===
    # Does M_σ's near-null cluster overlap with the plunge eigenspace?

    # The plunge eigenvectors in the whitened basis (where we computed eigenvalues)
    if len(plunge_indices) > 0:
        plunge_evecs_w = sorted_evecs[:, plunge_indices]
    else:
        plunge_evecs_w = np.zeros((N, 0))

    # Near-null eigenvectors in the whitened basis
    # (They are in sorted_Msig_evecs, which is from M_σ in G-whitened form)
    if near_null_count > 0:
        # Overlap: compute Gram matrix between near-null and plunge subspaces
        # In the whitened metric, this is just the SVD of their dot product
        if plunge_evecs_w.shape[1] > 0:
            # For each near-null eigenvector, project onto plunge eigenspace
            overlaps = np.zeros(near_null_count)
            for k in range(near_null_count):
                nk_vec = near_null_evecs[:, k]
                proj_norm = np.linalg.norm(plunge_evecs_w.T @ nk_vec)
                overlaps[k] = proj_norm
            max_overlap = overlaps.max() if len(overlaps) > 0 else 0.0
            mean_overlap = overlaps.mean() if len(overlaps) > 0 else 0.0
        else:
            overlaps = np.zeros(near_null_count)
            max_overlap = 0.0
            mean_overlap = 0.0
    else:
        overlaps = np.array([])
        max_overlap = 0.0
        mean_overlap = 0.0

    return {
        'N': N,
        'Ncap': Ncap,
        'L': L,
        'W': W,
        'prolate_spectrum': sorted_evals,
        'plunge_count': plunge_count,
        'plunge_evals': sorted_evals[plunge_indices] if len(plunge_indices) > 0 else np.array([]),
        'plunge_max_gap': plunge_max_gap,
        'lw_estimate': lw_estimate,
        'ev_Msig': sorted_Msig_evals,
        'near_null_count': near_null_count,
        'near_null_evals': near_null_evals,
        'Msig_max_gap': np.max(abs_gaps) if len(sorted_Msig_evals) > 1 else 0.0,
        'mu': mu,
        'nneg': nneg,
        'overlaps': overlaps,
        'max_overlap': max_overlap,
        'mean_overlap': mean_overlap,
    }

print("="*130)
print("ANGLE 2: Is M_σ's near-null cluster the Landau–Widom prolate plunge?")
print("="*130)

configs = [
    (14, 44, 4.0),
    (14, 44, 5.0),
    (14, 34, 4.0),
    (60, 90, 4.0),
    (120, 150, 4.0),
]

for (a, b, fourW) in configs:
    print(f"\n{'─'*130}")
    print(f"Configuration: [{a},{b}] 4W={fourW}")
    print(f"{'─'*130}")

    result = run(a, b, fourW)

    N = result['N']
    Ncap = result['Ncap']
    L = result['L']
    W = result['W']

    prolate_spec = result['prolate_spectrum']
    plunge_count = result['plunge_count']
    plunge_evals = result['plunge_evals']
    plunge_gap = result['plunge_max_gap']
    lw_est = result['lw_estimate']

    ev_Msig = result['ev_Msig']
    near_null_count = result['near_null_count']
    near_null_evals = result['near_null_evals']
    Msig_gap = result['Msig_max_gap']
    mu = result['mu']
    nneg = result['nneg']
    overlaps = result['overlaps']
    max_overlap = result['max_overlap']
    mean_overlap = result['mean_overlap']

    print(f"  Dimension: N={N} (Ncap={Ncap:.1f}), window length L={L}, bandwidth W={W}")

    # (1) Prolate concentration spectrum
    print(f"\n  [1] PROLATE CONCENTRATION SPECTRUM (generalized eigvals of C w.r.t. G):")
    print(f"      Eigenvalues (sorted descending, top 15):")
    top_k = min(15, len(prolate_spec))
    for i in range(top_k):
        print(f"        λ[{i:2d}] = {prolate_spec[i]:.6f}")

    print(f"      Largest spectral gap: {plunge_gap:.6f} (marks plunge start)")
    print(f"      Plunge region (after largest gap): {plunge_count} eigenvalues")
    if len(plunge_evals) > 0:
        print(f"      Plunge eigenvalues: {np.array2string(plunge_evals[:min(8,len(plunge_evals))], precision=4)}")

    print(f"      Landau–Widom plunge-width estimate: ln(2TB)/π² = ln({2*L*W:.1f})/π² ≈ {lw_est:.2f}")

    # (2) M_σ's near-null cluster
    print(f"\n  [2] M_σ's SPECTRUM (G-whitened):")
    print(f"      All eigenvalues (first 12): {np.array2string(ev_Msig[:min(12,len(ev_Msig))], precision=4)}")
    print(f"      Negative eigenvalues: {nneg}")
    if nneg > 0:
        print(f"        |λ_min| = {mu:.3e}")

    print(f"      Largest spectral gap in |λ|: {Msig_gap:.3e} (marks near-null cluster end)")
    print(f"      Near-null cluster (before gap): {near_null_count} eigenvalues")
    if len(near_null_evals) > 0:
        print(f"      Near-null eigenvalues: {np.array2string(near_null_evals[:min(8,len(near_null_evals))], precision=4)}")

    # (3) Comparison
    print(f"\n  [3] COMPARISON:")
    print(f"      Plunge width vs near-null count: {plunge_count} (plunge) vs {near_null_count} (near-null)")
    print(f"      Dimension match? (plunge ≈ near-null): {'LIKELY' if abs(plunge_count - near_null_count) <= 2 else 'UNLIKELY'}")

    # (4) Overlap with plunge eigenspace
    print(f"\n  [4] OVERLAP TEST (M_σ's negative eigenvec with plunge eigenspace):")
    if len(overlaps) > 0:
        print(f"      Max overlap (‖proj‖): {max_overlap:.4f}")
        print(f"      Mean overlap: {mean_overlap:.4f}")
        if max_overlap > 0.5:
            print(f"      ✓ STRONG overlap: M_σ's defect IS concentrated in the plunge")
        elif max_overlap > 0.2:
            print(f"      ~ MODERATE overlap: partial concentration in plunge")
        else:
            print(f"      ✗ WEAK overlap: defect not concentrated in plunge")
    else:
        print(f"      No plunge eigenvalues or no negative eigenvalue of M_σ")

print("\n" + "="*130)
print("SUMMARY:")
print("="*130)
print("""
HANDLE criteria (all must be true):
  1. Prolate spectrum shows clear cluster near 1, sharp plunge, cluster near 0
  2. Plunge width (count) ≈ Landau–Widom estimate ln(2TB)/π²
  3. Plunge width ≈ M_σ's near-null cluster dimension
  4. M_σ's negative eigenvector overlaps strongly (>0.5) with plunge eigenspace

If all four hold → the near-null cluster IS the prolate plunge (HANDLE: CONFIRMED)
If (1) & (2) hold but (3)-(4) fail → plunge exists but is NOT the M_σ defect (REPORT: plunge found, no match)
If (1)-(2) fail → the construction doesn't recover prolate spectrum (REPORT: negative result)
""")
print("="*130)
