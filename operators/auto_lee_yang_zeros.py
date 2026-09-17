import numpy as np
from mpmath import mp

def lee_yang_zeros_polynomial(n):
    # Placeholder for constructing a polynomial analogous to the Lee-Yang partition function
    coeffs = [mp.rand() for _ in range(n+1)]
    return np.poly1d(coeffs)

def main():
    mp.dps = 50
    n = 8  # Degree of the polynomial
    P = lee_yang_zeros_polynomial(n)
    roots = P.roots
    zeros_on_unit_circle = all(abs(root) == 1 for root in roots if abs(np.imag(root)) > 1e-5)
    print('ZEROS_ON_UNIT_CIRCLE' if zeros_on_unit_circle else 'ZEROS_NOT_ON_UNIT_CIRCLE')

if __name__ == '__main__':
    main()