"""S2878: plant off-line pair #a,#b (above node range) at Im = f*halfgap; report eps and SIGNED low-node errors."""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1]); dps=int(sys.argv[2]); N=int(sys.argv[3]); a,b=map(int,sys.argv[4].split(','))
W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True); S0=W02-WR-WP; M=S0.rows; mp.mp.dps=dps
nu=[2*mp.pi*n/L for n in range(-N,N+1)]
def term(z):
    w=4*mp.sin(z*L/2)**2/L; return mp.matrix([[w/((nu[i]+z)*(nu[j]+z)) for j in range(M)] for i in range(M)])
ga=mp.zetazero(a).imag; gb=mp.zetazero(b).imag; gm=(ga+gb)/2; h=(gb-ga)/2
base=term(ga)+term(-ga)+term(gb)+term(-gb)
G=[mp.zetazero(k).imag for k in range(1,5)]
for f in sys.argv[5:]:
    if f=='none': S=S0; lab='true zeta'
    else:
        z=gm+1j*mp.mpf(f)*h; T=term(z)+term(mp.conj(z))+term(-z)+term(-mp.conj(z))
        S=S0-base+mp.matrix([[mp.re(T[i,j]) for j in range(M)] for i in range(M)]); lab=f"Im={mp.nstr(mp.mpf(f)*h,3)}"
    E,V=mp.eigsy(S); o=sorted(range(M),key=lambda i:E[i]); xi=[V[r,o[0]] for r in range(M)]
    R=lambda zz: mp.fsum(xi[n]/(zz-nu[n]) for n in range(M))
    errs=[]
    for g in G:
        try: errs.append(mp.nstr(mp.findroot(R,g,tol=mp.mpf(10)**(-dps+10))-g,3))
        except Exception: errs.append('fail')
    print(f"{lab:>12}: eps={mp.nstr(E[o[0]],4):>11} #neg={sum(1 for e in E if e<0)} signed t_k-gamma_k(k=1..4)={errs}",flush=True)
