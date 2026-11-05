"""S2873: SIGNED node error t_k - gamma_k across N (one-sided convergence test)."""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1]); K=int(sys.argv[2])
for spec in sys.argv[3:]:
    N,dps=map(int,spec.split(':'))
    W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True); S=W02-WR-WP; M=S.rows; mp.mp.dps=dps
    nu=[2*mp.pi*n/L for n in range(-N,N+1)]
    E,V=mp.eigsy(S); i0=min(range(M),key=lambda i:E[i]); xi=[V[r,i0] for r in range(M)]
    R=lambda z: mp.fsum(xi[n]/(z-nu[n]) for n in range(M))
    out=[]
    for k in range(1,K+1):
        g=mp.zetazero(k).imag
        if g>nu[-1]*1.25: break
        t=mp.findroot(R,g,tol=mp.mpf(10)**(-dps+10)); out.append(mp.nstr(t-g,3))
    print(f"N={N} nu_N={mp.nstr(nu[-1],4)} signed t_k-gamma_k:",out,flush=True)
