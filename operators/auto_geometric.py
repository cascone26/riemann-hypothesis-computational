#!/usr/bin/env python3
import numpy as np
from scipy.linalg import eigvalsh
def geometric_structure_test(n):
    # Simulate a geometric structure matrix for the Riemann xi function
    G = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1):
            if i == j:
                G[i][j] = 2 * (-1)**i
            else:
                G[i][j] = G[j][i] = (-1)**(i+j) / (i - j)
    # Check if the geometric structure matrix is positive semi-definite
    min_eigenvalue = np.min(eigvalsh(G))
    return min_eigenvalue
if __name__ == '__main__':
    n = 10
    result = geometric_structure_test(n)
    print(f'Minimum Eigenvalue: {result}')