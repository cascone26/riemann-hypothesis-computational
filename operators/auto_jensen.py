import numpy as np
from scipy.special import gamma
def jensen_polynomial(n):
    # Construct Jensen polynomial coefficients.
    coeffs = [gamma(k + n / 2) / (gamma(k + 1) * gamma(n / 2)) for k in range(n + 1)]
    return np.poly1d(coeffs[::-1])
def jensen_positivity_test(degree):
    # Generate a Jensen polynomial and check its positivity properties.
    P = jensen_polynomial(degree)
    x_values = np.linspace(-1, 1, 1000)
    y_values = P(x_values)
    min_value = np.min(y_values)
    return min_value

degree = 8
result = jensen_positivity_test(degree)
print('Min value of Jensen polynomial:', result)