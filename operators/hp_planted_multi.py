"""S2867: planted off-line pairs, several at once / above node range. Args: lam2 dps N delta_frac pair1 pair2 ...
pair 'a,b' moves zeta zeros #a,#b to gm ± i*delta_frac*halfgap. Reports #neg, eigenvalues."""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1]); dps=int(sys.argv[2]); N=int(sys.argv[3]); df=mp.mpf(sys.argv[4]); pairs=[tuple(map(int,p.split(','))) for p in sys.argv[5:]]
W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True); S=W02-WR-WP; M=S.rows
nu=[2*mp.pi*n/L for n in range(-N,N+1)]
def term(z):
    w=4*mp.sin(z*L/2)**2/L; return mp.matrix([[w/((nu[i]+z)*(nu[j]+z)) for j in range(M)] for i in range(M)])
for (a,b) in pairs:
    ga=mp.zetazero(a).imag; gb=mp.zetazero(b).imag; gm=(ga+gb)/2; h=(gb-ga)/2; z=gm+1j*df*h
    T=term(z)+term(mp.conj(z))+term(-z)+term(-mp.conj(z))
    T=mp.matrix([[mp.re(T[i,j]) for j in range(M)] for i in range(M)])
    S=S-(term(ga)+term(-ga)+term(gb)+term(-gb))+T
    print(f"  planted pair #{a},#{b}: gm={mp.nstr(gm,6)} Im={mp.nstr(df*h,4)}")
E=sorted(mp.eigsy(S)[0])
print(f"lam2={lam2} N={N} nu_N={mp.nstr(nu[-1],5)}: #neg={sum(1 for e in E if e<0)}  lowest eigs={[mp.nstr(e,4) for e in E[:5]]}")
