#!/usr/bin/env python3
"""
§867++ PRINCIPAL SOS VECTOR ANALYSIS

Extract the principal SOS vector (the one associated with λ_min, which guarantees positivity
at the boundary). Analyze its structure, norms, and functional form across windows/bandwidths.

GOAL: Is there a pattern to the principal vector? Can we write it explicitly?
"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import cholesky, eigh
import warnings
warnings.filterwarnings('ignore')

def vonmangoldt(N):
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

NMAX = 3000
LAM = vonmangoldt(NMAX)
NS = np.arange(2, NMAX+1)
LOGN = np.log(NS)
LAMW = LAM[2:] / np.sqrt(NS)

def build_full_system(a, b, W, dr=0.004):
    """Build operators, Cholesky, and extract spectral info."""
    fourW = 4*W
    L_window = b - a
    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)

    Ncap = int((b-a)*2*W/np.pi)
    N = max(4, int(0.85*Ncap))
    CEN = np.linspace(a, b, N)

    def fejer(x):
        y = W*x
        out = np.ones_like(y)
        nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2
        return out

    Pm = np.array([fejer(RGRID-c) + fejer(RGRID+c) for c in CEN])
    G = (Pm @ Pm.T) * dr

    arch = (Pm * OMEGA[None, :]) @ Pm.T * dr / (2*np.pi)

    mask = LOGN < fourW
    ln = LOGN[mask]
    lw = LAMW[mask]
    cosu = np.cos(np.outer(ln, RGRID))
    prime = np.zeros((N, N))
    for j in range(N):
        gp = ((Pm[j][None, :]*Pm) @ cosu.T) * dr / (2*np.pi)
        prime[j, :] = 2.0 * (gp * lw[None, :]).sum(axis=1)
    prime = 0.5 * (prime + prime.T)

    def fc(c):
        y = W*(0.5j-c)
        y2 = W*(0.5j+c)
        return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2

    fp = np.array([fc(c) for c in CEN])
    pole = 2.0 * np.real(np.outer(fp, fp))

    Q = pole + arch - prime
    Q_sym = 0.5 * (Q + Q.T)

    # Eigendecomposition of Q w.r.t. G
    evs, vecs = eigh(Q_sym, G)

    # Cholesky
    try:
        L = cholesky(Q_sym, lower=True)
    except:
        return None

    # Principal SOS vector (associated with smallest eigenvalue)
    principal_ev = vecs[:, 0]  # Eigenvector for λ_min
    principal_l = L[0, :]      # First row of Cholesky factor

    return {
        'CEN': CEN,
        'Pm': Pm,
        'RGRID': RGRID,
        'OMEGA': OMEGA,
        'N': N,
        'margin': evs[0],
        'L': L,
        'principal_ev': principal_ev,
        'principal_l': principal_l,
        'evs': evs,
        'vecs': vecs,
        'G': G,
        'Q': Q_sym,
    }

def reconstruct_principal_vector_function(data, a, b, W):
    """Reconstruct the principal SOS vector as an explicit function."""
    CEN = data['CEN']
    Pm = data['Pm']
    RGRID = data['RGRID']
    principal_l = data['principal_l']

    # The principal SOS vector in Fejer basis: l_0 = sum_j L[0,j] * Fejer_j(r)
    # Approximate its values over [a,b]
    mask = (RGRID >= a) & (RGRID <= b)
    r_window = RGRID[mask]

    def fejer(x):
        y = W*x
        out = np.ones_like(y)
        nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2
        return out

    # Evaluate Fejer basis at window points
    Pm_eval = np.array([fejer(r_window-c) + fejer(r_window+c) for c in CEN])
    l0_over_window = principal_l @ Pm_eval

    return r_window, l0_over_window

def analyze_principal_vectors(test_cases):
    """Analyze principal vectors across multiple cases."""
    print("\n" + "="*90)
    print("PRINCIPAL SOS VECTOR ANALYSIS")
    print("="*90)

    results_by_window = {}

    for (a, b, fourW, label) in test_cases:
        W = fourW / 4
        data = build_full_system(a, b, W)

        if data is None:
            continue

        # Reconstruct the principal vector over the window
        r_win, l0_vals = reconstruct_principal_vector_function(data, a, b, W)

        # Store by window for comparison
        window_key = f"[{a},{b}]"
        if window_key not in results_by_window:
            results_by_window[window_key] = []

        results_by_window[window_key].append({
            'fourW': fourW,
            'label': label,
            'margin': data['margin'],
            'principal_l': data['principal_l'],
            'l0_norm': np.linalg.norm(data['principal_l']),
            'r_window': r_win,
            'l0_over_window': l0_vals,
            'CEN': data['CEN'],
            'N': data['N'],
        })

    # Print results by window
    for window, cases in sorted(results_by_window.items()):
        print(f"\n{window}:")
        print(f"  {'4W':>5} {'margin':>12} {'||l_0||':>10} {'max|l_0(r)|':>12} {'rms l_0':>12}")
        print("  " + "-"*65)

        for case in sorted(cases, key=lambda x: x['fourW']):
            l0_max = np.max(np.abs(case['l0_over_window']))
            l0_rms = np.sqrt(np.mean(case['l0_over_window']**2))
            print(f"  {case['fourW']:>5.1f} {case['margin']:>+12.4e} {case['l0_norm']:>10.4f} "
                  f"{l0_max:>12.4f} {l0_rms:>12.4f}")

            # Print the Fejer expansion of l_0
            print(f"    l_0 = ", end="")
            coeffs = case['principal_l']
            terms = []
            for j, c in enumerate(coeffs):
                if abs(c) > 1e-4:
                    sign = "+" if c > 0 else "-"
                    if j == 0:
                        terms.append(f"{c:+.4f}·φ_0")
                    else:
                        terms.append(f"{abs(c):+.4f}·φ_{j}")
            if terms:
                print(" ".join(terms[:3]))
                if len(terms) > 3:
                    print(f"    ... ({len(terms)} total terms)")
            else:
                print("(all coefficients < 1e-4)")

    # Cross-window comparison
    print("\n" + "="*90)
    print("CROSS-WINDOW STRUCTURE")
    print("="*90)

    for window, cases in sorted(results_by_window.items()):
        if len(cases) < 2:
            continue

        print(f"\n{window}: How does principal vector change with 4W?")
        for case in sorted(cases, key=lambda x: x['fourW']):
            l0_std = np.std(case['l0_over_window'])
            l0_mean = np.mean(case['l0_over_window'])
            print(f"  4W={case['fourW']:.1f}: mean(l_0)={l0_mean:+.4f}, std(l_0)={l0_std:.4f}, "
                  f"margin={case['margin']:+.4e}")

        # Test: is the principal vector shape stable?
        if len(cases) >= 2:
            l0_first = cases[0]['l0_over_window']
            l0_last = cases[-1]['l0_over_window']

            # Normalize and compare
            l0_first_norm = l0_first / (np.linalg.norm(l0_first) + 1e-15)
            l0_last_norm = l0_last / (np.linalg.norm(l0_last) + 1e-15)

            # Since they may be on different grids, just report their properties
            print(f"  Shape change: from μ={np.mean(l0_first):+.4f}, σ={np.std(l0_first):.4f} "
                  f"to μ={np.mean(l0_last):+.4f}, σ={np.std(l0_last):.4f}")

def main():
    print("="*90)
    print("§867++ PRINCIPAL SOS VECTOR EXTRACTION & STRUCTURAL ANALYSIS")
    print("="*90)

    test_cases = [
        (14, 44, 1.5, "[14,44], 4W=1.5"),
        (14, 44, 2.0, "[14,44], 4W=2.0"),
        (14, 44, 2.5, "[14,44], 4W=2.5"),
        (50, 80, 2.0, "[50,80], 4W=2.0"),
        (50, 80, 3.0, "[50,80], 4W=3.0"),
        (100, 130, 3.0, "[100,130], 4W=3.0"),
        (100, 130, 4.0, "[100,130], 4W=4.0"),
        (300, 330, 3.5, "[300,330], 4W=3.5"),
        (300, 330, 4.5, "[300,330], 4W=4.5"),
    ]

    analyze_principal_vectors(test_cases)

    print("\n" + "="*90)
    print("INTERPRETATION")
    print("="*90)
    print("""
1. PRINCIPAL VECTOR IDENTITY:
   The principal SOS vector l_0 (associated with λ_min) is the "most critical" one —
   it controls positivity at the boundary (the knife-edge). Its Fejer expansion
   l_0 = Σ_j L[0,j] φ_j(r) encodes the positivity structure geometrically.

2. STRUCTURAL STABILITY:
   If l_0 has a UNIVERSAL FORM (same structure at all heights/bandwidths after
   normalization), it would enable a unified reduction. The data above tests this.

3. FEJER COEFFICIENT PATTERNS:
   The magnitude of L[0,j] indicates which basis elements are important. Small
   coefficients (< 1e-4) can be dropped for a sparse representation.

4. NEXT: Build a SPARSE approximation of l_0 and test whether it generalizes
   across multiple (window, 4W) pairs. If it does, this is a step toward
   UNIFORM SOS.
""")

if __name__ == "__main__":
    main()
