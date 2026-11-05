"""S2874: Markov-type homotopy. tau = tau_K + tau_tail, tau_K = exact zero-side sum over first K=N zeros (+-),
tau_s = tau_K + s*tau_tail. At s=0 ground nodes are exactly +-gamma_1..gamma_N. Track t_k(s)-gamma_k."""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1]); dps=int(sys.argv[2]); N=int(sys.argv[3])
W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True); S=W02-WR-WP; M=S.rows; mp.mp.dps=dps
nu=[2*mp.pi*n/L for n in range(-N,N+1)]
G=[mp.zetazero(k).imag for k in range(1,N+1)]
def term(z):
    w=4*mp.sin(z*L/2)**2/L; return mp.matrix([[w/((nu[i]+z)*(nu[j]+z)) for j in range(M)] for i in range(M)])
TK=mp.matrix(M,M)
for g in G: TK+=term(g)+term(-g)
Tt=S-TK
print(f"lam2={lam2} N={N}: min eig tau_tail={mp.nstr(min(mp.eigsy(Tt)[0]),4)} (PSD expected under RH)")
for s in ['0','1e-30','1e-20','1e-10','1e-5','0.01','0.1','0.5','1']:
    Ts=TK+mp.mpf(s)*Tt
    E,V=mp.eigsy(Ts); i0=min(range(M),key=lambda i:E[i]); xi=[V[r,i0] for r in range(M)]
    R=lambda z: mp.fsum(xi[n]/(z-nu[n]) for n in range(M))
    out=[]
    for g in G[:6]:
        try: t=mp.findroot(R,g,tol=mp.mpf(10)**(-dps+10)); out.append(mp.nstr(t-g,3))
        except Exception: out.append('fail')
    print(f"  s={s:>6}: eps={mp.nstr(E[i0],4):>10}  t_k-gamma_k (k=1..6): {out}",flush=True)
