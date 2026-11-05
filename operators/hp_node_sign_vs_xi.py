"""Direct test: does the CCM/Shi ground vector agree in SIGN with Riemann's Xi at every node?
Full-index ground vector xi_j (j=-N..N) of tau=W02-WR-WP;  xihat(nu_k) = 2L^{-1/2} * (L/2) cos(pi k) xi_k
=> sign xihat(nu_k) = (-1)^k sign xi_k.  Compare with sign Xi(nu_k), Xi(t)=xi(1/2+it) (real)."""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1]); dps=int(sys.argv[2]); N=int(sys.argv[3])
W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True); S=W02-WR-WP; M=S.rows
E,V=mp.eigsy(S); i0=min(range(M),key=lambda i:E[i])
xi=[V[r,i0] for r in range(M)]
def Xi(t):
    s=mp.mpf(1)/2+1j*t
    return mp.re(s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s))
row=[]; agree=0
for k in range(0,N+1):
    a=(-1)**k*xi[N+k]; b=Xi(2*mp.pi*k/L)
    row.append('+' if a*b>0 else '-'); agree+= a*b>0
print(f"lam2={lam2} N={N}: sign(xihat(nu_k))*sign(Xi(nu_k)) for k=0..N: {''.join(row)}  agree={agree}/{N+1}")
# ratio profile: xihat(nu_k)/Xi(nu_k), normalized to k=0
r=[((-1)**k*xi[N+k])/Xi(2*mp.pi*k/L) for k in range(N+1)]
print("   ratio/ratio0 by k:", [mp.nstr(x/r[0],3) for x in r])
