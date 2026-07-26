#!/usr/bin/env python3
"""§314: PROVED worst-case Hermite margin (replaces §313B's sampled eps_d). Rigorous sufficient criterion:
He_d hyperbolic with simple roots => its values at its own critical points y_i (roots of He_{d-1}, since
He_d' = d He_{d-1}) alternate in sign. If a perturbation Dp satisfies |Dp(y_i)| < |He_d(y_i)| for all i,
then (He_d+Dp)(y_i) keeps those alternating signs at the FIXED points y_i, giving d sign changes on
(-inf,y_1),...,(y_{d-2},y_{d-1}),(y_{d-1},inf) => d real roots (leading coeff stays >0 for eps<1). QED-level.
For relative coefficient perturbation |Dp(y)| <= eps * sum_j |a_j||y|^j, so
   eps_d^proved = min_i |He_d(y_i)| / sum_j |a_j| |y_i|^j.
Compute for d=3..8, recompute N(d)=sqrt(d^3/(3 eps_d)), check coverage vs available direct checks (n<=96-d),
run any extra direct checks needed. Charisma: the bound = HANDLE (rigorous criterion, exact evaluation)."""
import mpmath as mp, pickle, os
mp.mp.dps=60
def hermite_coeffs(d):
    P=[mp.mpf(1)]; Q=[mp.mpf(0),mp.mpf(1)]
    if d==0: return P
    for k in range(1,d):
        R=[mp.mpf(0)]+Q
        for i in range(len(P)): R[i]-=k*P[i]
        P,Q=Q,R
    return Q  # low->high
def peval(co,x): return sum(c*x**j for j,c in enumerate(co))
print("PROVED margins eps_d = min_i |He_d(y_i)| / sum_j |a_j||y_i|^j  (y_i = roots of He_{d-1}):",flush=True)
epsp={}
for d in range(3,9):
    co=hermite_coeffs(d)
    cocrit=hermite_coeffs(d-1)
    ys=mp.polyroots(list(reversed(cocrit)),maxsteps=200,extraprec=100)
    ys=[mp.re(y) for y in ys]
    vals=[]
    for y in ys:
        num=abs(peval(co,y))
        den=sum(abs(c)*abs(y)**j for j,c in enumerate(co))
        vals.append(num/den)
    epsp[d]=min(vals)
    print(f"  d={d}: eps_proved = {mp.nstr(epsp[d],5)}   (sampled §313B was: {'1.0 .268 .127 .067 .038 .022'.split()[d-3]})",flush=True)
print("",flush=True)
print("re-derived thresholds N(d)=sqrt(d^3/(3 eps_proved)):",flush=True)
Nn={}
for d in range(3,9):
    Nn[d]=mp.sqrt(d**3/(3*epsp[d]))
    print(f"  d={d}: N_proved(d) ~ {mp.nstr(Nn[d],5)}",flush=True)
# coverage check vs 97 cached moments (J^{d,n} needs n+d <= 96)
gam=[mp.mpf(g) for g in pickle.load(open("operators/_gam96.pkl","rb"))]
def jensen_maximag(d,n):
    co=[mp.binomial(d,k)*gam[n+k] for k in range(d,-1,-1)]
    scale=mp.sqrt(abs(co[0])*abs(co[-1])); co=[c/scale for c in co]
    rts=mp.polyroots(co,maxsteps=300,extraprec=200)
    return max(abs(mp.im(r)) for r in rts)/max(abs(r) for r in rts)
print("",flush=True)
print("direct coverage to N_proved(d) using cached 97 moments:",flush=True)
for d in range(3,9):
    N=int(mp.ceil(Nn[d]))
    avail=96-d
    upto=min(N,avail)
    worst=mp.mpf(0)
    for n in range(0,upto+1):
        worst=max(worst,jensen_maximag(d,n))
    status="CLOSED for ALL n (proved-margin scheme)" if (N<=avail and worst<mp.mpf('1e-20')) else (f"checked n<={upto} ok; NEED moments to m~{N+d}" if worst<mp.mpf('1e-20') else "FAILURE?!")
    print(f"  d={d}: N={N}  checked n=0..{upto}  worst|Im|/|z|={mp.nstr(worst,3)}  => {status}",flush=True)
print("DONE",flush=True)
