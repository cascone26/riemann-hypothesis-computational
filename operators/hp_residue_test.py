"""S2873: direct test of F_N(gamma_k) * Xi'(gamma_k) * c / eps ~ const, c = F_N(0)/Xi(0)."""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1]); dps=int(sys.argv[2]); N=int(sys.argv[3]); K=int(sys.argv[4])
W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True); S=W02-WR-WP; M=S.rows; mp.mp.dps=dps
nu=[2*mp.pi*n/L for n in range(-N,N+1)]
E,V=mp.eigsy(S); i0=min(range(M),key=lambda i:E[i]); eps=E[i0]; xi=[V[r,i0] for r in range(M)]
F=lambda z: 2/mp.sqrt(L)*mp.sin(z*L/2)*mp.fsum(xi[n]/(z-nu[n]) for n in range(M))
def Xi(t):
    s=mp.mpf(1)/2+1j*t; return mp.re(s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s))
c=F(mp.mpf('1e-30'))/Xi(0)
print(f"lam2={lam2} N={N} eps={mp.nstr(eps,5)} c=F(0)/Xi(0)={mp.nstr(c,6)}")
for k in range(1,K+1):
    g=mp.zetazero(k).imag; d=mp.diff(Xi,g)
    print(f"  k={k} gamma={mp.nstr(g,8)}  F(gamma)={mp.nstr(F(g),4):>11}  r_k=F*Xi'*c/eps={mp.nstr(F(g)*d*c/eps,6)}")
