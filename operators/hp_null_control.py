"""Null control for S2863 sign-faithfulness: ground state of the Weil matrix WITHOUT primes (W02-WR) and
with primes SCALED (W02-WR-c*WP). Report sign pattern of xihat(nu_k) and agreement with Xi."""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1]); dps=int(sys.argv[2]); N=int(sys.argv[3])
W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True)
def Xi(t):
    s=mp.mpf(1)/2+1j*t; return mp.re(s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s))
X=[Xi(2*mp.pi*k/L) for k in range(N+1)]
for c in ['0','0.5','0.9','1','1.1']:
    S=W02-WR-mp.mpf(c)*WP; M=S.rows
    E,V=mp.eigsy(S); i0=min(range(M),key=lambda i:E[i]); xi=[V[r,i0] for r in range(M)]
    h=[(-1)**k*xi[N+k] for k in range(N+1)]
    sgn=''.join('+' if v>0 else '-' for v in h); agree=sum(1 for k in range(N+1) if h[k]*X[k]>0)
    sc=sum(1 for a,b in zip(h,h[1:]) if a*b<0)
    print(f"prime scale c={c}: e0={mp.nstr(E[i0],4)} xihat-sign@nodes={sgn} signchanges={sc} agree-with-Xi={max(agree,N+1-agree)}/{N+1}")
print("Xi sign @nodes:          ",''.join('+' if v>0 else '-' for v in X))
