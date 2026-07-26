#!/usr/bin/env python3
"""§315: the two-constant race — worst-case vs DIRECTIONAL Hermite margins.
The scheme's real perturbation is not adversarial: renormalization absorbs affine-in-k, leaving the cubic
direction g(k) = (k-d/2)^3 - (3(d/2)^2/4)(k-d/2)  (best-affine-reduced centered cubic), applied
multiplicatively coeff_k -> coeff_k * e^{t g(k)} with t = c_n/6 > 0 ONE-SIGNED (c_m>0 measured, §313A).
Measure:
 (1) proved worst-case eps_d (alternation criterion) extended to d=20 — decay law;
 (2) directional margins tau_d^+ and tau_d^- (max |t| keeping He_d-perturbed hyperbolic, per sign), d=3..20;
 (3) the race: d^3*tau_d vs d^3*eps_d-equivalents. If tau_d ~ const/d^3 (poly) while eps_d ~ q^d (geo), the
     effective ladder N(d)=sqrt(1/(3 tau_d)) is POLYNOMIAL — a structural improvement. Charisma: HANDLE (measured)."""
import mpmath as mp
mp.mp.dps=50
def hermite_coeffs(d):
    P=[mp.mpf(1)]; Q=[mp.mpf(0),mp.mpf(1)]
    if d==0: return P
    for k in range(1,d):
        R=[mp.mpf(0)]+Q
        for i in range(len(P)): R[i]-=k*P[i]
        P,Q=Q,R
    return Q  # low->high
def peval(co,x): return sum(c*x**j for j,c in enumerate(co))
def is_hyp(co_lowhigh):
    try:
        rts=mp.polyroots(list(reversed(co_lowhigh)),maxsteps=600,extraprec=300)
    except Exception:
        return False   # conditioning failure -> treat as not-hyperbolic (conservative for bisection)
    mxa=max(abs(r) for r in rts)
    return max(abs(mp.im(r)) for r in rts)/mxa < mp.mpf('1e-15')
print("(1) proved worst-case eps_d to d=20 (alternation bound), with ratio to previous:",flush=True)
prev=None
epsd={}
for d in range(3,21):
    co=hermite_coeffs(d); cocrit=hermite_coeffs(d-1)
    ys=[mp.re(y) for y in mp.polyroots(list(reversed(cocrit)),maxsteps=200,extraprec=100)]
    epsd[d]=min(abs(peval(co,y))/sum(abs(c)*abs(y)**j for j,c in enumerate(co)) for y in ys)
    print(f"  d={d:2d}: eps={mp.nstr(epsd[d],5):>10}  ratio={mp.nstr(epsd[d]/prev,4) if prev else '-'}",flush=True)
    prev=epsd[d]
print("",flush=True)
print("(2) directional cubic margins tau_d (+/-), g(k)=(k-d/2)^3-(3(d/2)^2/4)(k-d/2):",flush=True)
print(" d :   tau_d^+       tau_d^-      d^3*tau^+   d^3*tau^-   N_dir(d)=sqrt(1/(3 tau^+))",flush=True)
for d in range(3,21):
    co=hermite_coeffs(d)
    h=mp.mpf(d)/2
    def g(k): return (k-h)**3-(3*h*h/4)*(k-h)
    def hyp_at(t):
        co2=[c*mp.e**(t*g(k)) for k,c in enumerate(co)]
        return is_hyp(co2)
    taus={}
    for sgn in (1,-1):
        CAP=mp.mpf(64)
        lo,hi=mp.mpf(0),mp.mpf('0.5')
        while hyp_at(sgn*hi) and hi<CAP: hi*=2
        if hi>=CAP and hyp_at(sgn*CAP):
            taus[sgn]=CAP; continue
        for _ in range(50):
            mid=(lo+hi)/2
            if hyp_at(sgn*mid): lo=mid
            else: hi=mid
        taus[sgn]=lo
    Nd=mp.sqrt(1/(3*taus[1])) if taus[1]>0 else mp.inf
    print(f" {d:2d}: {mp.nstr(taus[1],5):>10} {mp.nstr(taus[-1],5):>10}  {mp.nstr(d**3*taus[1],5):>9} {mp.nstr(d**3*taus[-1],5):>9}   {mp.nstr(Nd,5)}",flush=True)
print("",flush=True)
print("READ: if d^3*tau stabilizes (poly law) while eps_d ratio stays ~0.6 (geometric), the race FLIPS:",flush=True)
print("      directional ladder N_dir(d) ~ d^{3/2}*const — polynomial, vs the exponential worst-case ladder.",flush=True)
print("DONE",flush=True)
