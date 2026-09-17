import numpy as np
from scipy.special import gamma
def jensen_polynomial_coefficients(d):
    # Coefficients of the d-th Jensen polynomial for the Riemann xi function
    return [(-1)**k * gamma(k + 0.5) / (gamma(0.5) * gamma(k + 1)) for k in range(d + 1)]
def main():
    d = 9  # Degree of the Jensen polynomial to test
    coeffs = jensen_polynomial_coefficients(d)
    P = np.poly1d(coeffs)
    roots = P.roots
    all_real_roots = all(np.isreal(root) for root in roots)
    print('ALL_REAL_ROOTS' if all_real_roots else 'NOT_ALL_REAL_ROOTS')
if __name__ == '__main__':
    main()