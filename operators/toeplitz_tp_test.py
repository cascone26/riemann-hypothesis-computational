#!/usr/bin/env python3
"""NEW METHOD (A14): total-positivity of the Xi/Jensen coefficient Toeplitz matrix.
RH <=> Xi in Laguerre-Polya <=> coeff sequence is Polya-Frequency (PF) <=> its Toeplitz
matrix [a_{i-j}] is Totally Positive (Aissen-Schoenberg-Whitney-Edrei). DIFFERENT from the
dead Hankel/moment test. Payoff: TP is closed under matrix products, and zeta = Euler product
=> if the Toeplitz matrix factors as a product of one-per-prime TP factors, TP (=RH) follows
from closure. Step 1: is the sequence even TP, and how far up the hierarchy?
Test via Neville elimination (Gasca-Pena): matrix is TP iff all Neville multipliers & pivots >=0."""
import mpmath as mp, json
mp.mp.dps=120
a=[mp.mpf(str(x)) for x in json.load(open('results/gorz_extended.json'))['gorz_a']]

def toeplitz_lower(seq,n):
    return [[seq[i-j] if i>=j else mp.mpf(0) for j in range(n)] for i in range(n)]

def neville_TP(M):
    """Neville elimination; returns (is_TP, first_bad, multipliers_min). TP iff all mults & pivots >=0."""
    n=len(M); A=[row[:] for row in M]; mults=[]; pivots=[]; bad=None
    for k in range(n):           # eliminate column k below diagonal using Neville (adjacent rows, bottom-up)
        for i in range(n-1,k,-1):
            if A[i-1][k]==0:
                if A[i][k]!=0: bad=('zero-pivot',i,k)
                continue
            m=A[i][k]/A[i-1][k]
            mults.append(float(m))
            if m<0 and bad is None: bad=('mult<0',i,k,float(m))
            for j in range(k,n):
                A[i][j]=A[i][j]-m*A[i-1][j]
    for k in range(n):
        pivots.append(float(A[k][k]))
        if A[k][k]<0 and bad is None: bad=('pivot<0',k,float(A[k][k]))
    is_TP = bad is None and all(p>=0 for p in pivots)
    return is_TP, bad, (min(mults) if mults else None), pivots

print("Xi/Jensen coeff sequence a_n (n=0..11):", [f"{float(x):.3e}" for x in a[:12]])
# TP2 (log-concavity) hierarchy first — clean, interpretable
print("\nTP2 ratios a_n^2/(a_{n-1}a_{n+1}) (>=1 => log-concave/TP2):")
r2=[float(a[n]**2/(a[n-1]*a[n+1])) for n in range(1,14)]
print("  ", [f"{x:.4f}" for x in r2])
print("  all TP2:", all(x>=1 for x in r2), " min ratio:", min(r2))

# full Neville TP test at growing sizes
print("\nNeville total-positivity test (lower-triangular Toeplitz):")
for n in [4,6,8,10,12,15]:
    if n>len(a): break
    M=toeplitz_lower(a,n)
    tp,bad,mmin,piv=neville_TP(M)
    print(f"  n={n:2d}: TP={tp}  min-multiplier={mmin:.4e}  first-violation={bad}")
