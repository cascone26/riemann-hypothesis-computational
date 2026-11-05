"""Planted off-line zero experiment (S2865). Exact explicit-side tau (dps high) + zero-side surgery:
replace on-line zeros {±ga, ±gb} by {±(gm+r), ±(gm-r)}, r = sqrt(sigma) (sigma<0 => off-line pair gm ± i|r|).
term(z)_mn = 4 sin^2(zL/2) / (L (nu_m+z)(nu_n+z)).  Track e0 and node-sign agreement of xihat with
the planted Xi'(t) = Xi(t) * (t^2-z1^2)(t^2-z2^2)/((t^2-ga^2)(t^2-gb^2))."""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1]); dps=int(sys.argv[2]); N=int(sys.argv[3]); ia=int(sys.argv[4]); ib=int(sys.argv[5])
W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True); S0=W02-WR-WP; M=S0.rows
idx=list(range(-N,N+1)); nu=[2*mp.pi*n/L for n in idx]
ga=mp.zetazero(ia).imag; gb=mp.zetazero(ib).imag; gm=(ga+gb)/2; h=(gb-ga)/2
def term(z):
    T=mp.matrix(M,M); w=4*mp.sin(z*L/2)**2/L
    for i in range(M):
        for j in range(M): T[i,j]=w/((nu[i]+z)*(nu[j]+z))
    return T
def Xi(t):
    s=mp.mpf(1)/2+1j*t; return mp.re(s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s))
base=term(ga)+term(-ga)+term(gb)+term(-gb)
Xn=[Xi(2*mp.pi*k/L) for k in range(N+1)]
print(f"lam2={lam2} N={N} planting on zeros #{ia},#{ib}: ga={mp.nstr(ga,8)} gb={mp.nstr(gb,8)} gm={mp.nstr(gm,8)} half-gap={mp.nstr(h,5)}")
for frac in ['1','0.5','0.1','0','-0.01','-0.1','-0.3','-1']:
    sig=mp.mpf(frac)*h**2; r=mp.sqrt(sig) if sig>=0 else 1j*mp.sqrt(-sig)
    zs=[gm+r,gm-r,-(gm+r),-(gm-r)]
    new=mp.matrix(M,M)
    for z in zs: new+=term(z)
    new=mp.matrix([[mp.re(new[i,j]) for j in range(M)] for i in range(M)])
    S=S0-base+new
    E,V=mp.eigsy(S); order=sorted(range(M),key=lambda i:E[i]); i0=order[0]
    xi=[V[r_,i0] for r_ in range(M)]
    hh=[(-1)**k*xi[N+k] for k in range(N+1)]
    Xp=[Xn[k]*mp.re(((nu[N+k]**2-zs[0]**2)*(nu[N+k]**2-zs[1]**2))/((nu[N+k]**2-ga**2)*(nu[N+k]**2-gb**2))) for k in range(N+1)]
    ag=sum(1 for k in range(N+1) if hh[k]*Xp[k]>0); ag=max(ag,N+1-ag)
    print(f" sigma/h^2={frac:>5}: e0={mp.nstr(E[i0],5):>12} e1={mp.nstr(E[order[1]],5):>12} #neg={sum(1 for e in E if e<0)}  sign-agree w/ planted Xi'={ag}/{N+1}  xihat signs={''.join('+' if v>0 else '-' for v in hh)}")
