#!/usr/bin/env python3
"""Extract and characterize the NEGATIVE subspace of the indefinite archimedean+boundary form A+B.

§11 found: A+B = W02-WR is INDEFINITE (n_neg negative eigenvalues, growing with c), and tau=(A+B)-Wp is
PD only because the prime term Wp neutralizes that negative subspace. So the negative directions of A+B are
the seat of the whole problem. This pulls them out and asks:
  1. What ARE they? frequency profile (n -> omega=2pi n/L), localized or spread?
  2. Do they lock onto zeta zeros (unlike the full near-null mode, which did not)?
  3. How does the prime form act there? sign & size of v^T Wp v on each negative mode.
     (For tau=(A+B)-Wp to stay PD on v, need v^T(A+B)v - v^T Wp v >= 0, i.e. v^T Wp v <= v^T(A+B)v < 0:
      the prime term must be EVEN MORE negative on these modes. We measure whether/how that holds.)

A+B's negative eigenvalues are O(1)-scale (not the tiny tau margin), so double-precision eigh on the
rigorously-assembled (Arb) matrix midpoints is sufficient for characterization.
"""
import sys, os, math, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "external", "groskin_2607.02828", "code"))
import numpy as np
from flint import arb, ctx
from arb_ldlt_certify import arb_closed_forms, arb_kappa, arb_J, prime_powers_up_to

def build_split_np(c, N, prec=2000):
    """Return (arch, prime) as float64 numpy arrays (rigorously assembled in Arb, midpoints taken)."""
    ctx.prec = prec
    S, CC, XC, L = arb_closed_forms(N, c, prec)
    PI = arb.pi(); sp2 = 16*PI*PI; l2 = L*L
    pref02 = 32*L*(L/4).sinh()**2; kappa = arb_kappa(L); J = arb_J(L)
    pdata = prime_powers_up_to(c)
    weights = [arb(p).log()*(arb(q)**arb('-0.5')) for (q,p) in pdata]
    positions = [arb(q).log() for (q,p) in pdata]; M=len(weights)
    def Ss(nn): return S[nn] if nn>=0 else -S[-nn]
    DIM=2*N+1; Aa=np.zeros((DIM,DIM)); Ap=np.zeros((DIM,DIM))
    for i in range(DIM):
        n=i-N
        for j in range(i,DIM):
            m=j-N
            num=l2-sp2*m*n; den=(l2+sp2*m*m)*(l2+sp2*n*n); W02=pref02*num/den
            if n==m: WR=kappa+2*CC[abs(n)]+J-(2/L)*XC[abs(n)]
            else:    WR=(Ss(m)-Ss(n))/(PI*(n-m))
            Wp=arb(0)
            for idx in range(M):
                y=positions[idx]
                if n==m: q=2*(1-y/L)*(2*PI*n*y/L).cos()
                else:    q=((2*PI*m*y/L).sin()-(2*PI*n*y/L).sin())/(PI*(n-m))
                Wp+=weights[idx]*q
            a=float((W02-WR).mid()); p=float(Wp.mid())
            Aa[i,j]=a; Aa[j,i]=a; Ap[i,j]=p; Ap[j,i]=p
    return Aa, Ap, float(L.mid())

ZEROS=[14.1347,21.0220,25.0109,30.4249,32.9351,37.5862,40.9187,43.3271,48.0052,49.7738]

def profile(v, N, L):
    mass=[(i-N, v[i]*v[i]) for i in range(len(v))]
    mass.sort(key=lambda t:-t[1]); tot=sum(m for _,m in mass) or 1
    top=[(n, m/tot, 2*math.pi*n/L) for n,m in mass[:5]]
    # dominant |omega| and nearest zero
    domn = mass[0][0]; domomega = 2*math.pi*domn/L
    znear=min(ZEROS,key=lambda g:abs(g-abs(domomega)))
    return top, domomega, znear, sum(m for _,m in mass[:5])/tot

def main():
    N=64
    print(f"=== negative subspace of A+B (archimedean+boundary), N={N} ===")
    out=[]
    for c in [14, 50, 200]:
        Aa,Ap,L=build_split_np(c,N)
        w,V=np.linalg.eigh(Aa)              # ascending eigenvalues
        neg=[k for k in range(len(w)) if w[k]<0]
        print(f"\nc={c:4d} (L={L:.3f}): A+B has {len(neg)} negative eigenvalues; "
              f"most-negative {w[neg[0]] if neg else float('nan'):+.4e}")
        modes=[]
        for k in neg:
            v=V[:,k]; lam=w[k]
            top,domomega,znear,top5=profile(v,N,L)
            vApv=float(v@Ap@v)              # prime form on this arch-negative mode
            vAav=float(v@Aa@v)              # = lam (check)
            taum=vAav-vApv                  # tau value on v (not the min, just on this mode)
            near = abs(abs(domomega)-znear)<0.7
            print(f"  eig#{k}: lam(A+B)={lam:+.4e}  domfreq |w|={abs(domomega):6.3f}"
                  f"{'  ~ZERO '+str(znear) if near else '':<14}  top5mass={100*top5:.0f}%")
            print(f"          v^T Wp v = {vApv:+.4e}   v^T tau v = {taum:+.4e}  "
                  f"({'prime MORE negative (restores PD)' if vApv<lam else 'prime does NOT cover' if taum<0 else 'ok'})")
            modes.append(dict(k=k,lam_AB=lam,dom_omega=domomega,nearest_zero=znear,locks_zero=bool(near),
                              top5mass=top5,vWpv=vApv,tau_on_v=taum))
        out.append(dict(c=c,L=L,n_neg=len(neg),modes=modes))
    json.dump(out,open("results/neg_subspace_AB.json","w"),indent=2)
    # summary: do negative modes lock onto zeros? does prime always cover them?
    print("\n=== READ ===")
    total_modes=sum(len(r['modes']) for r in out)
    locked=sum(1 for r in out for m in r['modes'] if m['locks_zero'])
    covered=sum(1 for r in out for m in r['modes'] if m['tau_on_v']>=0)
    print(f"  negative A+B modes examined: {total_modes}")
    print(f"  lock onto a zeta zero (|domfreq - gamma|<0.7): {locked}/{total_modes}")
    print(f"  prime term keeps tau>=0 on the mode: {covered}/{total_modes}")
    print("wrote results/neg_subspace_AB.json")

if __name__=="__main__":
    main()
