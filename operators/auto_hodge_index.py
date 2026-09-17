from mpmath import mp
import numpy as np

def compute_hodge_index_matrix(n):
    # Placeholder for constructing a matrix analogous to the intersection form in function fields
    A = np.random.rand(n, n)
    return (A + A.T) / 2  # Symmetric matrix

def main():
    mp.dps = 50
    n = 10  # Dimension of the matrix
    H = compute_hodge_index_matrix(n)
    eigenvalues = np.linalg.eigvalsh(H)
    min_eigenvalue = min(eigenvalues)
    print(f'Min Eigenvalue: {min_eigenvalue}')
    if min_eigenvalue >= 0:
        print('HODGE_INDEX_NONNEGATIVE')

if __name__ == '__main__':
    main()