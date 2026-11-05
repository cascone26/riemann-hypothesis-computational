"""Envelope m_N(nu_k)=xihat_N(nu_k)/Xi(nu_k) across N (fixed lambda): local convergence test of CCM xihat_N -> c Xi.
Also report eigenvalue e0, min envelope, and sign agreement. Saves JSON."""
import sys, json, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1]); out={}
def Xi(t):
    s=mp.mpf(1)/2+1j*t
    return mp.re(s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s))
for spec in sys.argv[2:]:
    N,dps=map(int,spec.split(':'))
    W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True); S=W02-WR-WP; M=S.rows
    E,V=mp.eigsy(S); i0=min(range(M),key=lambda i:E[i]); xi=[V[r,i0] for r in range(M)]
    nu=[2*mp.pi*k/L for k in range(N+1)]
    r=[((-1)**k*xi[N+k])/Xi(nu[k]) for k in range(N+1)]
    m=[x/r[0] for x in r]
    out[N]={'e0':mp.nstr(E[i0],6),'nu':[float(v) for v in nu],'m':[float(x) for x in m]}
    print(f"N={N} dps={dps} e0={mp.nstr(E[i0],5)} sign-agree={sum(1 for x in m if x>0)}/{N+1} min m={mp.nstr(min(m),4)} at k={m.index(min(m))}",flush=True)
    print("   m:",[mp.nstr(x,3) for x in m],flush=True)
json.dump(out,open(f'results/hp_center/envelope_l{int(lam2)}.json','w'))
