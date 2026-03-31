"""
Verify our LMFDB data against mpmath's independent zero computation.
This confirms our pipeline is correct before we start analysis.
"""

import mpmath
import os

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def load_zeros(path=None, count=None):
    """Load zeros from our data file as mpmath mpf values."""
    if path is None:
        # Prefer LMFDB high-precision data, fall back to Odlyzko
        lmfdb = os.path.join(DATA_DIR, "zeros_100k.txt")
        odlyzko = os.path.join(DATA_DIR, "zeros_odlyzko_100k")
        path = lmfdb if os.path.exists(lmfdb) else odlyzko

    zeros = []
    with open(path, "r") as f:
        for i, line in enumerate(f):
            if count and i >= count:
                break
            line = line.strip()
            if line:
                zeros.append(mpmath.mpf(line))
    return zeros


def verify_against_mpmath(n_verify=100):
    """
    Verify first n zeros against mpmath.zetazero().
    mpmath computes these independently using the Euler-Maclaurin formula.
    """
    mpmath.mp.dps = 35  # 35 decimal places (our data has 31)

    our_zeros = load_zeros(count=n_verify)
    if len(our_zeros) < n_verify:
        print(f"Only have {len(our_zeros)} zeros loaded. Verifying those.")
        n_verify = len(our_zeros)

    print(f"Verifying first {n_verify} zeros against mpmath...")
    print(f"Working precision: {mpmath.mp.dps} decimal places")
    print("-" * 70)

    max_error = mpmath.mpf(0)
    errors = []

    for i in range(1, n_verify + 1):
        mp_zero = mpmath.zetazero(i)
        mp_imag = mp_zero.imag
        our_val = our_zeros[i - 1]
        error = abs(mp_imag - our_val)
        errors.append(float(error))

        if error > max_error:
            max_error = error

        if i <= 10 or i % 20 == 0:
            print(f"Zero #{i:>4d}: LMFDB = {mpmath.nstr(our_val, 20)} | "
                  f"mpmath = {mpmath.nstr(mp_imag, 20)} | "
                  f"error = {mpmath.nstr(error, 5)}")

    print("-" * 70)
    print(f"Max absolute error: {mpmath.nstr(max_error, 10)}")
    print(f"Mean absolute error: {mpmath.nstr(sum(errors) / len(errors), 10)}")

    if max_error < mpmath.mpf("1e-25"):
        print("\nVERIFICATION PASSED — data matches to 25+ decimal places.")
        return True
    else:
        print(f"\nWARNING — max error {max_error} exceeds threshold.")
        return False


if __name__ == "__main__":
    verify_against_mpmath(100)
