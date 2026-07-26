"""§255 prong: continuous-omega envelope sweep (YellowBot's proposal, room 2026-07-24).
Same psi as _h_decay.py (c=13 exact-Q near-null). Instead of sampling only at prime
powers, sweep omega densely in (0,1] with UNIT weight w=1: t=(e^{(1-om)L},(1-om)L,1).
Outputs log10|h_env(om)|, the sin-prefactor-corrected log10|h/om| (YellowBot's flag),
and fits: power (log om), exp(-kappa/om) (1/om), Gaussian (om^2), plus local slopes.
"""
import sys, os, mpmath as mp
VENV=os.path.expanduser("~/projects/riemann/cvs_venv/lib/python3.14/site-packages")
sys.path.insert(0, VENV); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import connes_cvs.operator as O
from cvs_fast import precompute_components
c=9; N=40; dps=120; T=300
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
print(f"c={c} L={mp.nstr(L,5)} lam_min log10={float(mp.log10(abs(eigs[k0]))):.2f}", flush=True)
# continuous sweep: omega in (0,1], unit weight. Prime omegas for c=13 span ~[0.03,0.74];
# go a bit past both ends to see the envelope shape cleanly.
M=25
oms=[mp.mpf('0.001')*mp.power(mp.mpf(50),mp.mpf(i)/mp.mpf(24)) for i in range(25)]
print(f"{'omega':>8} {'log10|h_env|':>13} {'log10|h/om|':>12}", flush=True)
rows=[]
for om in oms:
    logn=(1-om)*L; nn=mp.e**logn
    t=(nn, logn, mp.mpf(1))
    Vv={m:O.psi_prime(mp.mpf(m),L,[t]) for m in idx}
    Vd={m:O.psi_prime_deriv(mp.mpf(m),L,[t]) for m in idx}
    h=(psi.T*assemble(Vv,Vd)*psi)[0]
    omf=float(om); lh=float(mp.log10(abs(h))) if h!=0 else float('-inf')
    lc=lh-float(mp.log10(om))
    rows.append((omf,lh,lc))
    print(f"{omf:8.4f} {lh:13.3f} {lc:12.3f}", flush=True)
import math
def fit(x,y):
    n_=len(x); mx=sum(x)/n_; my=sum(y)/n_
    b=sum((xi-mx)*(yi-my) for xi,yi in zip(x,y))/sum((xi-mx)**2 for xi in x)
    r=sum((xi-mx)*(yi-my) for xi,yi in zip(x,y))/math.sqrt(sum((xi-mx)**2 for xi in x)*sum((yi-my)**2 for yi in y))
    return b, my-b*mx, r
# fit the small-omega tail (lower half), prefactor-corrected column
tail=[r for r in rows if r[0]<0.4 and r[1]>-90]
x_log=[math.log10(r[0]) for r in tail]; x_inv=[1/r[0] for r in tail]; x_sq=[r[0]**2 for r in tail]
y=[r[2] for r in tail]
for name,xs in (("POWER  y~b*log10(om)",x_log),("EXPINV y~b*(1/om)",x_inv),("GAUSS  y~b*om^2",x_sq)):
    b,a,r=fit(xs,y); print(f"{name}: b={b:.3f} a={a:.2f} r={r:.5f}", flush=True)
# local slopes d log10|h/om| / d log10 om, and vs 1/om (PurpleBot's regression, corrected)
print("\nlocal slopes (prefactor-corrected):")
for i in range(1,len(rows)):
    o0,_,c0=rows[i-1]; o1,_,c1=rows[i]
    if c0>-90 and c1>-90 and o1<0.5:
        s=(c1-c0)/(math.log10(o1)-math.log10(o0))
        print(f"  om~{(o0+o1)/2:.3f} slope={s:.1f} 1/om={2/(o0+o1):.1f}", flush=True)
