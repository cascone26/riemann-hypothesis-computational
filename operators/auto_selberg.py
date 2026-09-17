from scipy.sparse import diags
from numpy.linalg import eigh
import numpy as np

def selberg_positivity_test(n):
    # Construct a sparse matrix mimicking the Selberg zeta Laplacian spectrum.
    data = [np.ones(n-1), -2*np.ones(n), np.ones(n-1)]
    offsets = [-1, 0, 1]
    L = diags(data, offsets).toarray()
    P = (L + L.T) / 2
    A = np.eye(n)
    B = np.random.rand(n, n)
    B = (B + B.T) / 2
    tau = A + B - P
    eigenvalues, _ = eigh(tau)
    return min(eigenvalues)

result = selberg_positivity_test(10)
print('Min eigenvalue of tau:', result)