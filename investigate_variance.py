"""
Investigate the spacing variance discrepancy: we got 0.16, GUE predicts 0.286.
This script checks whether it's a normalization issue or a real signal.
"""

import math
import os

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def load_zeros(count=None):
    path = os.path.join(DATA_DIR, "zeros_odlyzko_100k")
    zeros = []
    with open(path, "r") as f:
        for i, line in enumerate(f):
            if count and i >= count:
                break
            line = line.strip()
            if line:
                zeros.append(float(line))
    return zeros


def test_normalization_methods():
    """Try different normalization approaches and see which gives GUE variance."""
    zeros = load_zeros()
    print(f"Loaded {len(zeros)} zeros\n")

    # Raw spacings
    raw = [zeros[i + 1] - zeros[i] for i in range(len(zeros) - 1)]

    # Method 1: Our original — local density normalization
    method1 = []
    for i in range(len(zeros) - 1):
        t = zeros[i]
        density = math.log(t / (2 * math.pi)) / (2 * math.pi)
        method1.append(raw[i] * density)

    # Method 2: Global mean normalization (divide by mean spacing)
    mean_raw = sum(raw) / len(raw)
    method2 = [s / mean_raw for s in raw]

    # Method 3: Local mean normalization (window of 100 neighbors)
    method3 = []
    window = 50
    for i in range(len(raw)):
        lo = max(0, i - window)
        hi = min(len(raw), i + window)
        local_mean = sum(raw[lo:hi]) / (hi - lo)
        method3.append(raw[i] / local_mean)

    # Method 4: Using the unfolding formula directly
    # N(t) ~ (t/2pi) * log(t/2pi) - t/2pi + 7/8
    # Unfolded zeros: u_n = N(t_n), then spacings of u_n should be ~1
    def N_approx(t):
        return (t / (2 * math.pi)) * math.log(t / (2 * math.pi)) - t / (2 * math.pi) + 7 / 8

    unfolded = [N_approx(t) for t in zeros]
    method4 = [unfolded[i + 1] - unfolded[i] for i in range(len(unfolded) - 1)]

    for name, spacings in [("Method 1: local density", method1),
                           ("Method 2: global mean", method2),
                           ("Method 3: local window", method3),
                           ("Method 4: smooth unfolding", method4)]:
        n = len(spacings)
        mean = sum(spacings) / n
        var = sum((s - mean) ** 2 for s in spacings) / n
        print(f"{name}")
        print(f"  Mean: {mean:.8f}  Variance: {var:.8f}")
        print(f"  (GUE expects mean~1.0, variance~0.286)")
        print()


def convergence_test():
    """Check if variance converges to GUE as we use more zeros."""
    zeros = load_zeros()

    print("Convergence test (Method 4: smooth unfolding)")
    print("-" * 50)

    def N_approx(t):
        return (t / (2 * math.pi)) * math.log(t / (2 * math.pi)) - t / (2 * math.pi) + 7 / 8

    for n_zeros in [100, 500, 1000, 5000, 10000, 50000, 100000]:
        subset = zeros[:n_zeros]
        unfolded = [N_approx(t) for t in subset]
        spacings = [unfolded[i + 1] - unfolded[i] for i in range(len(unfolded) - 1)]

        mean = sum(spacings) / len(spacings)
        var = sum((s - mean) ** 2 for s in spacings) / len(spacings)
        print(f"  N={n_zeros:>6d}: mean={mean:.6f}, var={var:.6f}")


if __name__ == "__main__":
    print("=" * 60)
    print("VARIANCE DISCREPANCY INVESTIGATION")
    print("=" * 60)
    print()
    test_normalization_methods()
    print()
    convergence_test()
