"""Oscillation-structure probe beyond the TP threshold: sign-change counts of every eigenvector of each
parity block of tau (= Shi S_NL).  Oscillation matrix => eigenvector of the j-th LARGEST eigenvalue has
exactly j-1 sign changes (ground state the most oscillatory)."""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
from hp_tp_test import sectors
lam2=float(sys.argv[1]); dps=int(sys.argv[2])
def sc(v,tol):
    s=[x for x in v if abs(x)>tol]; return sum(1 for a,b in zip(s,s[1:]) if a*b<0)
for N in [int(a) for a in sys.argv[3:]]:
    W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True); S=W02-WR-WP
    Ev,Od=sectors(S,N)
    for name,B in [('even',Ev),('odd',Od)]:
        E,V=mp.eigsy(B); n=B.rows; order=sorted(range(n),key=lambda i:-E[i])  # largest first
        counts=[sc([V[r,i] for r in range(n)],mp.mpf(10)**(-dps//3)) for i in order]
        print(f"N={N} {name}: dim={n} sign changes (largest->smallest eig): {counts}  ideal={list(range(n))}")
