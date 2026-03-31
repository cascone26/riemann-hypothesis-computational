"""
Sanity check: is our WKB "winner" just re-discovering the smooth
counting function, or does it actually capture oscillatory corrections?

If N_WKB(t_n) ≈ N_smooth(t_n) for all n, the search found nothing new.
If N_WKB(t_n) tracks N_exact(t_n) = n better than N_smooth does, there's
genuine content beyond the smooth part.
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import math
import os

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
PLOTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")


def load_zeros(count=200):
    path = os.path.join(DATA_DIR, "zeros_100k.txt")
    zeros = []
    with open(path, "r") as f:
        for i, line in enumerate(f):
            if i >= count:
                break
            line = line.strip()
            if line:
                zeros.append(float(line))
    return np.array(zeros)


ALPHA, BETA, A, B = 0.876416, -0.167183, 4.628555, 19.160354

def N_smooth(E):
    if E <= 2:
        return 0.0
    return (E / (2*math.pi)) * math.log(E / (2*math.pi)) - E / (2*math.pi) + 7/8

def N_wkb(E):
    def integrand(x):
        U = A * x**ALPHA
        V = B * x**BETA
        disc = E**2 - 4*U*V
        if disc <= 0 or U <= 0:
            return 0.0
        return math.sqrt(disc) / U
    try:
        result, _ = quad(integrand, 0.01, 200, limit=200, epsabs=1e-8)
        return result / (2 * math.pi)
    except:
        return 0.0


def main():
    zeros = load_zeros(100)

    print("SANITY CHECK: Is our WKB winner just the smooth function?")
    print("=" * 65)

    print(f"\n{'n':>3} {'zero':>10} {'N_exact':>8} {'N_smooth':>10} {'N_WKB':>10} "
          f"{'|sm-ex|':>8} {'|wkb-ex|':>8} {'WKB better?':>12}")
    print("-" * 75)

    smooth_errors = []
    wkb_errors = []
    wkb_wins = 0

    for i in range(min(50, len(zeros))):
        E = zeros[i]
        n_exact = i + 1
        n_smooth = N_smooth(E)
        n_wkb = N_wkb(E)

        err_smooth = abs(n_smooth - n_exact)
        err_wkb = abs(n_wkb - n_exact)
        smooth_errors.append(err_smooth)
        wkb_errors.append(err_wkb)

        better = "YES" if err_wkb < err_smooth else "no"
        if err_wkb < err_smooth:
            wkb_wins += 1

        if i < 30 or i % 10 == 0:
            print(f"{n_exact:>3} {E:>10.4f} {n_exact:>8d} {n_smooth:>10.4f} {n_wkb:>10.4f} "
                  f"{err_smooth:>8.4f} {err_wkb:>8.4f} {better:>12}")

    print(f"\n{'='*65}")
    print(f"Mean |N_smooth - n|:  {np.mean(smooth_errors):.6f}")
    print(f"Mean |N_WKB - n|:     {np.mean(wkb_errors):.6f}")
    print(f"WKB wins:             {wkb_wins}/{len(smooth_errors)}")
    print(f"Ratio (WKB/smooth):   {np.mean(wkb_errors)/np.mean(smooth_errors):.4f}")

    if np.mean(wkb_errors) < np.mean(smooth_errors) * 0.8:
        print("\n>>> WKB is GENUINELY better than smooth (>20% improvement)")
        print(">>> The search found something beyond the smooth counting function.")
    elif np.mean(wkb_errors) < np.mean(smooth_errors):
        print("\n>>> WKB is slightly better but within noise.")
        print(">>> Probably just a better fit to the smooth function, not oscillatory content.")
    else:
        print("\n>>> WKB is NOT better than smooth.")
        print(">>> The search just re-discovered the smooth counting function.")

    # Plot the residuals N(t_n) - n for both
    fig, ax = plt.subplots(figsize=(14, 6))
    ns = range(1, len(smooth_errors) + 1)
    smooth_resid = [N_smooth(zeros[i]) - (i+1) for i in range(len(smooth_errors))]
    wkb_resid = [N_wkb(zeros[i]) - (i+1) for i in range(len(wkb_errors))]

    ax.plot(ns, smooth_resid, 'r.-', markersize=3, label=f'N_smooth - n (MAE={np.mean(smooth_errors):.3f})')
    ax.plot(ns, wkb_resid, 'b.-', markersize=3, label=f'N_WKB - n (MAE={np.mean(wkb_errors):.3f})')
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.set_xlabel('n')
    ax.set_ylabel('N(t_n) - n')
    ax.set_title('Residuals: Smooth vs WKB Counting Functions at Zeta Zeros')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'sanity_check.png'), dpi=150)
    plt.close()
    print(f"\nPlot saved.")


if __name__ == "__main__":
    main()
