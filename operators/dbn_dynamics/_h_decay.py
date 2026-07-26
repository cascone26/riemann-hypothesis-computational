import sys, os, mpmath as mp
VENV=os.path.expanduser("~/projects/riemann/cvs_venv/lib/python3.14/site-packages")
sys.path.insert(0, VENV); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import connes_cvs.operator as O
from cvs_fast import precompute_components
c=13; N=40; dps=120; T=300
mp.mp.dps=dps
comp=precompute_components(c,N,T,dps); idx=comp["idx"]; L=comp["L"]; DIM=2*N+1
Vp=mp.matrix(DIM,N+1); Vp[N,0]=mp.mpf(1); r2=1/mp.sqrt(2)
for k in range(1,N+1): Vp[N+k,k]=r2; Vp[N-k,k]=r2
arc={n:comp["arcpol"][n] for n in idx}; arcd={n:comp["arcpol_d"][n] for n in idx}
def assemble(vm,dm):
    Q=mp.matrix(DIM,DIM)
    for i in range(DIM):
        m=i-N
        for j in range(DIM):
            nn=j-N; Q[i,j]=dm[nn] if m==nn else (vm[m]-vm[nn])/(m-nn)
    for i in range(DIM):
        for j in range(i+1,DIM):
            a=(Q[i,j]+Q[j,i])/2; Q[i,j]=a; Q[j,i]=a
    return Vp.T*Q*Vp
pdata,_=O.prime_powers_up_to(c)
fv={}; fd={}
for n in idx:
    fv[n]=arc[n]+sum(O.psi_prime(mp.mpf(n),L,[t]) for t in pdata)
    fd[n]=arcd[n]+sum(O.psi_prime_deriv(mp.mpf(n),L,[t]) for t in pdata)
eigs,vecs=mp.eigsy(assemble(fv,fd)); k0=min(range(N+1),key=lambda i:eigs[i])
psi=mp.matrix([vecs[r,k0] for r in range(N+1)])
print(f"c={c} L={mp.nstr(L,5)} lam_min log10={float(mp.log10(abs(eigs[k0]))):.2f}")
print(f"{'n':>3} {'omega_n':>9} {'log10|h|':>10} {'log10(w=L/sqrtn)':>16} {'log10|overlap|':>15}")
rows=[]
for (n,logn,w) in pdata:
    Vv={nn:O.psi_prime(mp.mpf(nn),L,[(n,logn,w)]) for nn in idx}
    Vd={nn:O.psi_prime_deriv(mp.mpf(nn),L,[(n,logn,w)]) for nn in idx}
    h=(psi.T*assemble(Vv,Vd)*psi)[0]
    om=float(1-logn/L); lh=float(mp.log10(abs(h))); lw=float(mp.log10(w))
    rows.append((n,om,lh,lw,lh-lw))
    print(f"{n:>3} {om:9.4f} {lh:10.2f} {lw:16.2f} {lh-lw:15.2f}")
# fit log|h| vs log(omega) and vs omega
import math
om=[r[1] for r in rows if r[1]>0]; lh=[r[2] for r in rows if r[1]>0]
n_=len(om); 
def fit(x,y):
    mx=sum(x)/n_; my=sum(y)/n_
    b=sum((xi-mx)*(yi-my) for xi,yi in zip(x,y))/sum((xi-mx)**2 for xi in x)
    return b, my-b*mx
b1,a1=fit([math.log(o) for o in om], lh)   # power law log|h| = a + b log omega
b2,a2=fit(om, lh)                            # exp law log|h| = a + b omega
print(f"\nPOWER: log10|h| ~ {b1:.2f}*log10(omega) + {a1:.2f}  => h ~ omega^{b1:.1f}  (b~c? c={c})")
print(f"EXP:   log10|h| ~ {b2:.2f}*omega + {a2:.2f}")
