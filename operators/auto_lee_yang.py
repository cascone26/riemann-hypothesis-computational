import numpy as np
from mpmath import *
mp.dps = 50

def lee_yang_positivity_test(n):
    # Generate a matrix similar to the Lee-Yang zeros distribution.
    angles = [2 * pi * k / n for k in range(1, n + 1)]
    roots = [exp(j*angle) for angle in angles]
    A = np.diag([sum(root**k for root in roots) for k in range(n)])
    B = eye(n)
    P = (A + A.T) / 2
    tau = A + B - P
    return min(eig(tau))

result = lee_yang_positivity_test(10)
print('Min eigenvalue of tau:', result)