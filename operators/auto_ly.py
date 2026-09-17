#!/usr/bin/env python3
import numpy as np
def lee
lee_yang_test(N=100):
    # Create a matrix representing interactions
   
 J = np.random.rand(N, N)
    H = np.dot(J, J.T)  # Ensure the matrix is p
positive semi-definite
    eigvals = np.linalg.eigvalsh(H)
    min_eigval
min_eigval = np.min(eigvals)
    print(min_eigval)
    return min_eigval 
>= 0
if __name__ == '__main__':
    lee_yang_test()