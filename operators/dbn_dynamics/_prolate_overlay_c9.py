"""GreenBot, §255 kill-shot scaffold (room 2026-07-24): PARAMETER-FREE arch-mode overlay.
Question (YellowBot's kinematics-vs-arithmetic): is the near-null psi's envelope pure
arch/prolate kinematics? Test WITHOUT fits: diagonalize the ARCH-ONLY operator Qa
(the prolate side per LEDGER (SS)252), expand the full near-null psi in Qa's eigenbasis,
and (a) print the coefficient spectrum |<phi_k|psi>| vs arch eigenvalue mu_k,
(b) reconstruct the continuous envelope h_env(om) from psi truncated to its top-K arch
modes at TWO depths K (PurpleBot's truncation null: only where both K agree with each
other but disagree with full psi is there structure beyond truncation).
If psi ~= a single deep arch mode + tiny corrections -> envelope is arch kinematics,
ALL arithmetic lives in the ratio==1 placement (the pre-registered branch).
Config matches (SS)253b: c=9 N=30 dps=100 (cheap, ~2 min class).
"""
import sys, os, mpmath as mp
VENV=os.path.expanduser("~/projects/riemann/cvs_venv/lib/python3.14/site-packages")
sys.path.insert(0, VENV); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import connes_cvs.operator as O
from cvs_fast import precompute_components
c=9; N=30; dps=100; T=250
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
# full operator -> near-null psi
pdata,_=O.prime_powers_up_to(c)
fv={}; fd={}
for n in idx:
    fv[n]=arc[n]+sum(O.psi_prime(mp.mpf(n),L,[t]) for t in pdata)
    fd[n]=arcd[n]+sum(O.psi_prime_deriv(mp.mpf(n),L,[t]) for t in pdata)
eigs,vecs=mp.eigsy(assemble(fv,fd)); k0=min(range(N+1),key=lambda i:eigs[i])
psi=mp.matrix([vecs[r,k0] for r in range(N+1)])
print(f"c={c} L={mp.nstr(L,5)} lam_min log10={float(mp.log10(abs(eigs[k0]))):.2f}", flush=True)
# arch-only diagonalization (the prolate side)
Qa=assemble(arc,arcd)
aeigs,avecs=mp.eigsy(Qa)
order=sorted(range(N+1),key=lambda i:float(aeigs[i]))
print("\n(a) ARCH-MODE DECOMPOSITION of psi  (k sorted by arch eigenvalue mu_k)", flush=True)
print(f"{'k':>3} {'log10|mu_k|':>12} {'sgn':>4} {'log10|<phi_k|psi>|':>19}")
coeffs=[]
for rk,i in enumerate(order):
    phi=mp.matrix([avecs[r,i] for r in range(N+1)])
    ck=(phi.T*psi)[0]; coeffs.append((rk,i,ck))
    lmu=float(mp.log10(abs(aeigs[i]))) if aeigs[i]!=0 else float('-inf')
    lck=float(mp.log10(abs(ck))) if ck!=0 else float('-inf')
    sg='+' if aeigs[i]>0 else '-'
    print(f"{rk:>3} {lmu:12.3f} {sg:>4} {lck:19.3f}", flush=True)
# (b) two-truncation envelope reconstruction
def env_at(psi_v, om):
    logn=(1-om)*L; nn=mp.e**logn
    t=(nn, logn, mp.mpf(1))
    Vv={m:O.psi_prime(mp.mpf(m),L,[t]) for m in idx}
    Vd={m:O.psi_prime_deriv(mp.mpf(m),L,[t]) for m in idx}
    return (psi_v.T*assemble(Vv,Vd)*psi_v)[0]
def truncate(K):
    byc=sorted(coeffs,key=lambda x:-abs(x[2]))[:K]
    v=mp.matrix(N+1,1)
    for _,i,ck in byc:
        phi=mp.matrix([avecs[r,i] for r in range(N+1)])
        for r in range(N+1): v[r]+=ck*phi[r]
    nrm=mp.sqrt((v.T*v)[0])
    for r in range(N+1): v[r]/=nrm
    return v
psiK1=truncate(8); psiK2=truncate(16)
oms=[mp.mpf(x)/mp.mpf(100) for x in (2,4,7,12,20,30,45,60,75,90)]
print("\n(b) ENVELOPE: full psi vs top-8 vs top-16 arch-mode truncations  (log10|h_env|)", flush=True)
print(f"{'omega':>6} {'full':>9} {'K=8':>9} {'K=16':>9}")
for om in oms:
    hf=env_at(psi,om); h1=env_at(psiK1,om); h2=env_at(psiK2,om)
    f=lambda h: f"{float(mp.log10(abs(h))):9.2f}" if h!=0 else "     -inf"
    print(f"{float(om):6.2f} {f(hf)} {f(h1)} {f(h2)}", flush=True)
print("\nREAD: coeff spectrum ~single dominant mode -> kinematics; K=8~K=16~full at all om ->", flush=True)
print("envelope is low-arch-mode content (kinematics); K's agree but differ from full -> beyond-arch structure.", flush=True)
