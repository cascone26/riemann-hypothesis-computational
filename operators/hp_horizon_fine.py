"""S2879: detection horizon Im_crit(N) for an off-line pair planted at zeros #a,#b (bisection on lambda_min<0)."""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1]); a,b=map(int,sys.argv[2].split(','))
for spec in sys.argv[3:]:
    N,dps=map(int,spec.split(':'))
    W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True); S0=W02-WR-WP; M=S0.rows; mp.mp.dps=dps
    nu=[2*mp.pi*n/L for n in range(-N,N+1)]
    def term(z):
        w=4*mp.sin(z*L/2)**2/L; return mp.matrix([[w/((nu[i]+z)*(nu[j]+z)) for j in range(M)] for i in range(M)])
    ga=mp.zetazero(a).imag; gb=mp.zetazero(b).imag; gm=(ga+gb)/2
    base=S0-(term(ga)+term(-ga)+term(gb)+term(-gb))
    def lmin(y):
        z=gm+1j*y; T=term(z)+term(mp.conj(z))+term(-z)+term(-mp.conj(z))
        S=base+mp.matrix([[mp.re(T[i,j]) for j in range(M)] for i in range(M)])
        return min(mp.eigsy(S)[0])
    lo,hi=mp.mpf('1e-15'),mp.mpf('0.01')
    if lmin(hi)>=0: print(f"N={N}: not detected even at Im=3",flush=True); continue
    if lmin(lo)<0: print(f"N={N}: detected already at Im=0.001",flush=True); continue
    for it in range(60):
        mid=mp.sqrt(lo*hi)
        if lmin(mid)<0: hi=mid
        else: lo=mid
    print(f"lam2={lam2} N={N} nu_N={mp.nstr(nu[-1],5)} gamma_m={mp.nstr(gm,6)}: Im_crit ~ {mp.nstr((lo+hi)/2,5)}",flush=True)
