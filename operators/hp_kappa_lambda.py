"""S2876: universality of kappa = |t_k-g_k|*|Xi'(g_k)|^2/eps across lambda."""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
def Xi(t):
    s=mp.mpf(1)/2+1j*t; return mp.re(s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s))
for spec in sys.argv[1:]:
    lam2,N,dps=spec.split(':'); lam2=float(lam2); N=int(N); dps=int(dps)
    W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True); S=W02-WR-WP; M=S.rows; mp.mp.dps=dps
    nu=[2*mp.pi*n/L for n in range(-N,N+1)]
    E,V=mp.eigsy(S); i0=min(range(M),key=lambda i:E[i]); eps=E[i0]; xi=[V[r,i0] for r in range(M)]
    R=lambda z: mp.fsum(xi[n]/(z-nu[n]) for n in range(M))
    ks=[]
    for k in range(1,6):
        g=mp.zetazero(k).imag
        if g>nu[-1]: break
        t=mp.findroot(R,g,tol=mp.mpf(10)**(-dps+10)); d=mp.diff(Xi,g)
        ks.append(mp.nstr((t-g)*d**2/eps,4))
    print(f"lam2={lam2} N={N} eps={mp.nstr(eps,4)}: kappa_k={ks}",flush=True)
