"""S2871: node error |t_k-gamma_k| for k=1..6 vs N (rate fit for T1)."""
import sys, json, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1]); res={}
def polymul(a,b):
    r=[mp.mpf(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]+=x*y
    return r
for spec in sys.argv[2:]:
    N,dps=map(int,spec.split(':'))
    W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True); S=W02-WR-WP; M=S.rows
    mp.mp.dps=dps; G=[mp.zetazero(k).imag for k in range(1,11)]
    nu=[2*mp.pi*n/L for n in range(-N,N+1)]
    E,V=mp.eigsy(S); i0=min(range(M),key=lambda i:E[i]); xi=[V[r,i0] for r in range(M)]
    # Newton-refine zeros of R near each gamma (avoid full polyroots)
    R=lambda z: mp.fsum(xi[n]/(z-nu[n]) for n in range(M))
    errs=[]
    for g in G:
        try: t=mp.findroot(R,g,tol=mp.mpf(10)**(-dps+10)); errs.append(abs(t-g))
        except Exception as e: errs.append(mp.nan)
    res[N]={'eps':mp.nstr(E[i0],5),'err':[mp.nstr(x,4) for x in errs]}
    print(f"N={N} dps={dps} nu_N={mp.nstr(nu[-1],5)} eps={mp.nstr(E[i0],4)} errs:",[mp.nstr(x,3) for x in errs],flush=True)
json.dump(res,open('results/hp_center/rate_scan_l14.json','w'))
