#!/usr/bin/env python3
"""VECTOR A: (i) does the diagonal self-energy floor min_ω τ(ω) stay bounded below as c→∞ (necessary
condition for the uniform inequality Q_A ≥ Q_P), and (ii) the large-ω asymptotic of the archimedean
self-energy Q_A(ω) (needed for any analytic estimate). Uses the validated continuous closed forms.
"""
import sys, os, math, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "external", "groskin_2607.02828", "code"))
import mpmath as mp
from arb_ldlt_certify import prime_powers_up_to

def funcs(c, dps=40):
    mp.mp.dps=dps; L=mp.log(c); PI=mp.pi
    pdata=prime_powers_up_to(c)
    W=[(mp.log(p)/mp.sqrt(q), mp.log(q)) for (q,p) in pdata]
    psi_q=mp.digamma(mp.mpf(1)/4)
    def geom(n):
        w=2*PI*n/L; w2=w*w; gS=gCC=gX1=gX2=mp.mpf(0); k=0
        while True:
            ck=2*k+mp.mpf('0.5'); e=mp.e**(-ck*L); den=ck*ck+w2
            gS+=e/den
            if n!=0: gCC+=e*w2/(ck*den)
            gX1+=e*ck/den; gX2+=e*(ck*ck-w2)/(den*den)
            if e<mp.mpf(10)**(-dps-5) and k>2: break
            k+=1
        return gS,gCC,gX1,gX2
    eL=mp.e**L; U=L/2
    J=-2*mp.log(mp.e**U+1)+mp.log(mp.e**(2*U)+1)+2*mp.atan(mp.e**U)+mp.log(2)-PI/2
    kappa=mp.log(4*PI*(eL-1)/(eL+1))+mp.euler
    pref02=32*L*mp.sinh(L/4)**2; sp2=16*PI*PI; l2=L*L
    def QA(omega):
        n=omega*L/(2*PI); z=mp.mpf(1)/4+1j*PI*n/L
        psi=mp.digamma(z); psi1=mp.polygamma(1,z)
        gS,gCC,gX1,gX2=geom(n)
        CC=-(mp.re(psi)-psi_q)/2+gCC; XC=mp.re(psi1)/4-L*gX1-gX2
        WR=kappa+2*CC+J-(2/L)*XC
        num=l2-sp2*n*n; den=(l2+sp2*n*n)**2; W02=pref02*num/den
        return float(W02-WR)
    def QP(omega):
        s=mp.mpf(0)
        for (wt,y) in W: s+=wt*2*(1-y/L)*mp.cos(omega*y)
        return float(s)
    return QA,QP,float(L)

def main():
    print("=== A: diagonal floor min_ω τ(ω) vs c, and Q_A large-ω asymptotic ===")
    print(f"{'c':>6} {'L':>7} {'min_ω τ (floor)':>16} {'argmin ω':>10} {'Q_A(40)':>10} {'Q_A(80)':>10}")
    rows=[]
    for c in [50,200,800,3200,12800]:
        QA,QP,L=funcs(c)
        # scan floor over a broad band avoiding the zero-resonance spikes (take between-zero minima)
        best=None
        w=0.5
        while w<=40:
            t=QA(w)-QP(w)
            if best is None or t<best[1]: best=(w,t)
            w+=0.05
        rows.append(dict(c=c,L=L,floor=best[1],argmin=best[0],QA40=QA(40),QA80=QA(80)))
        print(f"{c:>6} {L:>7.3f} {best[1]:>16.5f} {best[0]:>10.3f} {QA(40):>10.3f} {QA(80):>10.3f}")
    # Q_A asymptotic: fit Q_A(ω) ~ a*log(ω)+b on a big grid at largest c
    QA,QP,L=funcs(12800)
    ws=[float(x) for x in range(20,120,5)]
    ys=[QA(w) for w in ws]
    import statistics
    lx=[math.log(w) for w in ws]
    n=len(lx);sx=sum(lx);sy=sum(ys);sxx=sum(x*x for x in lx);sxy=sum(x*y for x,y in zip(lx,ys))
    a=(n*sxy-sx*sy)/(n*sxx-sx*sx); b=(sy-a*sx)/n
    print(f"\n  Q_A(ω) ~ {a:.4f}·ln(ω) + {b:.4f}  (archimedean self-energy grows logarithmically)")
    json.dump({"floor_vs_c":rows,"QA_asymptotic":{"a_lnomega":a,"b":b}},
              open("results/A_uniform_floor.json","w"),indent=2)
    print("\n=== VERDICT (A) ===")
    fl=[r['floor'] for r in rows]
    print(f"  diagonal floor across c: {[f'{x:.4f}' for x in fl]}")
    if all(x>0 for x in fl):
        trend = fl[-1]/fl[0]
        print(f"  floor stays POSITIVE; ratio last/first = {trend:.2f} "
              f"({'~bounded below (uniform inequality plausible)' if trend>0.3 else 'shrinking toward 0 — inequality gets tight'})")
    else:
        print("  floor goes NEGATIVE somewhere — diagonal inequality already fails (would be a surprise).")
    print("wrote results/A_uniform_floor.json")

if __name__=="__main__":
    main()
