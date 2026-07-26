"""§255: direct measurement of the odd-moment suppression ladder of the near-null.
E(omega) = -(1/pi) sum_k (-1)^k (2 pi omega)^{2k+1}/(2k+1)! * M_k,
M_k = psi^T Q[x^{2k+1}] psi (same divided-difference assembly as psi_prime).
Deep-tail sweeps show E(om) ~ c0*om with c0 ~ 1e-34 at c=9: the ladder is
suppressed, not killed. Print log10|M_k| + sign for k=0..24, then self-check by
resumming the series at two omegas against the direct assembly.
Usage: python _moment_ladder.py <c> <N> <dps>
"""
import sys, os, mpmath as mp
VENV=os.path.expanduser("~/projects/riemann/cvs_venv/lib/python3.14/site-packages")
sys.path.insert(0, VENV); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import connes_cvs.operator as O
from cvs_fast import precompute_components
c=int(sys.argv[1]); N=int(sys.argv[2]); dps=int(sys.argv[3]); T=300
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
print(f"c={c} N={N} dps={dps} L={mp.nstr(L,5)} lam_min log10={float(mp.log10(abs(eigs[k0]))):.2f}", flush=True)
KMAX=25
Ms=[]
print(f"{'k':>3} {'2k+1':>5} {'log10|M_k|':>11} {'sign':>5} {'log10|series coeff at om=1|':>26}", flush=True)
for k in range(KMAX):
    p=2*k+1
    vm={m: mp.mpf(m)**p for m in idx}
    dm={m: p*mp.mpf(m)**(p-1) for m in idx}
    M=(psi.T*assemble(vm,dm)*psi)[0]
    Ms.append(M)
    coef=abs(M)*(2*mp.pi)**p/mp.factorial(p)/mp.pi
    print(f"{k:>3} {p:>5} {float(mp.log10(abs(M))):11.3f} {'+' if M>0 else '-':>5} {float(mp.log10(coef)):26.3f}", flush=True)
# self-check: resum series vs direct sin assembly at two omegas
for om in (mp.mpf('0.001'), mp.mpf('0.01'), mp.mpf('0.05')):
    S=-sum((-1)**k*(2*mp.pi*om)**(2*k+1)/mp.factorial(2*k+1)*Ms[k] for k in range(KMAX))/mp.pi
    t=(mp.e**((1-om)*L),(1-om)*L,mp.mpf(1))
    Vv={m:O.psi_prime(mp.mpf(m),L,[t]) for m in idx}
    Vd={m:O.psi_prime_deriv(mp.mpf(m),L,[t]) for m in idx}
    h=(psi.T*assemble(Vv,Vd)*psi)[0]
    print(f"om={float(om):.3f}: series log10={float(mp.log10(abs(S))):.3f} direct log10={float(mp.log10(abs(h))):.3f} rel.err={float(abs((S-h)/h)):.2e}", flush=True)
