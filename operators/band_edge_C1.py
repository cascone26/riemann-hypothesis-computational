#!/usr/bin/env python3
"""RUNG 3 / VECTOR C1 (make-or-break): is the negative-eigenvalue band FIXED in frequency as c->inf?

§12 hint: A+B's negative modes are the low harmonics; their max frequency was ~4.8,6.4,5.9 at c=14,50,200
(roughly constant) while n_neg=4,7,10 grew. If the band EDGE omega*(c)=max negative-mode frequency
converges to a constant while n_neg~const*L grows only from mode density, then:
  - the negative part is a FIXED compact frequency band (core), the rest is tail,
  - the "finite core + controlled tail" reduction (RH_SOLVE_PLAN §1) is VIABLE.
If omega* -> infinity, the reduction as stated fails and we need a different frame.

Also captures A1: the diagonal self-energies Q_A(omega)=arch[i,i], Q_P(omega)=Wp[i,i] vs omega=2pi n/L.
Double-precision eigh is fine (arch's negative eigenvalues are O(1); we only need signs + which frequency).
"""
import sys, os, math, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "external", "groskin_2607.02828", "code"))
import numpy as np
from flint import arb, ctx
from arb_ldlt_certify import arb_closed_forms, arb_kappa, arb_J, prime_powers_up_to

def build(c, N, prec=1200):
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

def dom_freq(v, N, L):
    k=int(np.argmax(v*v)); return 2*math.pi*(k-N)/L

def main():
    N=80
    print(f"=== C1: negative-band edge omega*(c) as c grows (N={N}) ===")
    print(f"{'c':>6} {'L=lnc':>7} {'n_neg':>6} {'n_neg/L':>8} {'omega* (band edge)':>18} {'min eig':>12}")
    rows=[]
    for c in [14, 50, 200, 800, 3200, 12800]:
        Aa,Ap,L=build(c,N)
        w,V=np.linalg.eigh(Aa)
        neg=[k for k in range(len(w)) if w[k]<0]
        freqs=[abs(dom_freq(V[:,k],N,L)) for k in neg]
        omega_star=max(freqs) if freqs else 0.0
        rows.append(dict(c=c,L=L,n_neg=len(neg),omega_star=omega_star,min_eig=float(w[0])))
        print(f"{c:>6} {L:>7.3f} {len(neg):>6} {len(neg)/L:>8.3f} {omega_star:>18.4f} {w[0]:>12.3e}")
    # A1: diagonal self-energies vs frequency, at the largest c
    Aa,Ap,L=build(12800,N)
    diag=[]
    for i in range(2*N+1):
        n=i-N; diag.append((2*math.pi*n/L, Aa[i,i], Ap[i,i], Aa[i,i]-Ap[i,i]))
    json.dump(dict(band=rows, diag_self_energy_c12800=diag),
              open("results/band_edge_C1.json","w"), indent=2)
    # verdict
    print("\n=== VERDICT (C1) ===")
    os_ = [r['omega_star'] for r in rows]
    nl = [r['n_neg']/r['L'] for r in rows]
    conv_band = (max(os_[-3:]) - min(os_[-3:])) < 1.5   # last 3 within 1.5 in omega
    flat_dens = (max(nl[-3:]) - min(nl[-3:])) < 0.5     # n_neg/L roughly constant
    print(f"  omega* last-3 spread = {max(os_[-3:])-min(os_[-3:]):.2f}  (converged if <1.5): {conv_band}")
    print(f"  n_neg/L last-3 spread = {max(nl[-3:])-min(nl[-3:]):.3f} (flat if <0.5): {flat_dens}")
    if conv_band and flat_dens:
        print("  => PASS: negative part is a FIXED frequency band; n_neg grows ONLY from mode density (~L).")
        print("     The 'finite core + controlled tail' reduction is VIABLE. Rung 3 supported.")
    elif conv_band:
        print("  => PARTIAL: band edge converges but density law unclear; reduction likely viable, verify.")
    else:
        print("  => FAIL as stated: band edge grows with c; need a different frame (rescale/structure fn).")
    print("wrote results/band_edge_C1.json")

if __name__=="__main__":
    main()
