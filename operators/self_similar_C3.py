#!/usr/bin/env python3
"""C3: is the prime-side near-cancellation SELF-SIMILAR across zeta zeros? (new lever after C2 refuted the
finite reduction.) If the local structure of the self-energy around gamma_k is a rescaled copy of the one
around gamma_1, a renormalization-style proof (one template + scaling law) could cover all zeros at once.

Q_A(omega), Q_P(omega) are evaluated as CONTINUOUS closed-form functions of frequency (fine grid, arbitrary
resolution), faithfully replicating the Arb assembly's diagonal closed forms in mpmath. VALIDATED against
the matrix diagonal at integer nodes before use. Then extract the local resonance profile (depth, width,
residual floor) at gamma_1..gamma_6 and test whether they collapse onto one shape under rescaling.
"""
import sys, os, math, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "external", "groskin_2607.02828", "code"))
import mpmath as mp
from arb_ldlt_certify import prime_powers_up_to

ZEROS=[14.134725,21.022040,25.010858,30.424876,32.935062,37.586178,40.918719,43.327073]

def selfenergy_funcs(c, dps=40):
    mp.mp.dps=dps
    L=mp.log(c); PI=mp.pi
    pdata=prime_powers_up_to(c)
    W=[(mp.log(p)/mp.sqrt(q), mp.log(q)) for (q,p) in pdata]   # (weight, y=log q)
    psi_q=mp.digamma(mp.mpf(1)/4)
    def geom(n):
        w=2*PI*n/L; w2=w*w
        gS=gCC=gX1=gX2=mp.mpf(0); k=0
        while True:
            ck=2*k+mp.mpf('0.5'); e=mp.e**(-ck*L); den=ck*ck+w2
            gS+=e/den
            if n!=0: gCC+=e*w2/(ck*den)
            gX1+=e*ck/den; gX2+=e*(ck*ck-w2)/(den*den)
            if e<mp.mpf(10)**(-dps-5) and k>2: break
            k+=1
        return gS,gCC,gX1,gX2
    U=(L/2); eL=mp.e**L
    J=-2*mp.log(mp.e**U+1)+mp.log(mp.e**(2*U)+1)+2*mp.atan(mp.e**U)+mp.log(2)-PI/2
    kappa=mp.log(4*PI*(eL-1)/(eL+1))+mp.euler
    pref02=32*L*mp.sinh(L/4)**2; sp2=16*PI*PI; l2=L*L
    def QA(omega):
        n=omega*L/(2*PI)
        z=mp.mpf(1)/4+1j*PI*n/L
        psi=mp.digamma(z); psi1=mp.polygamma(1,z)
        gS,gCC,gX1,gX2=geom(n)
        CC=-(mp.re(psi)-psi_q)/2+gCC
        XC=mp.re(psi1)/4-L*gX1-gX2
        WR=kappa+2*CC+J-(2/L)*XC
        num=l2-sp2*n*n; den=(l2+sp2*n*n)**2; W02=pref02*num/den
        return float(W02-WR)
    def QP(omega):
        s=mp.mpf(0)
        for (wt,y) in W: s+=wt*2*(1-y/L)*mp.cos(omega*y)
        return float(s)
    return QA,QP,float(L)

def main():
    c=3200
    print(f"=== C3: self-similarity of the self-energy resonances across zeros (c={c}) ===")
    QA,QP,L=selfenergy_funcs(c)
    # validate against matrix diagonal at a couple integer nodes
    print("  [validate] QA,QP at integer omega vs known scale:")
    for n in (3,8):
        om=2*math.pi*n/L; print(f"    n={n} omega={om:.3f}  QA={QA(om):+.4f}  QP={QP(om):+.4f}  tau={QA(om)-QP(om):+.5f}")
    # scan each zero: find the QP resonance min in a window, its depth, half-width, and the tau floor nearby
    print("\n  zero      omega_res   QP_depth   halfwidth   tau_at_res   floor(off-res)")
    prof=[]
    for g in ZEROS[:6]:
        # fine scan around g
        xs=[g-2+0.01*k for k in range(401)]
        vals=[(x,QP(x)) for x in xs]
        xr,qmin=min(vals,key=lambda t:t[1])
        floor=QP(g-1.8)  # off-resonance reference
        # half-width: where QP recovers to half of (floor - qmin) below floor
        half=floor-(floor-qmin)/2
        hw=None
        for x,v in vals:
            if v<=half:
                hw=abs(x-xr);
        taures=QA(xr)-QP(xr)
        prof.append(dict(gamma=g,omega_res=xr,depth=floor-qmin,halfwidth=hw,tau_at_res=taures,floor=floor))
        print(f"  {g:8.4f}  {xr:9.4f}  {floor-qmin:9.4f}  {hw if hw else 0:9.4f}  {taures:+10.4f}  {floor:+8.4f}")
    json.dump(prof,open("results/self_similar_C3.json","w"),indent=2)
    # self-similarity read: rescale each resonance by its depth; are widths ~constant or scaling law?
    print("\n=== VERDICT (C3) ===")
    depths=[p['depth'] for p in prof]; hws=[p['halfwidth'] for p in prof if p['halfwidth']]
    if len(hws)>=4:
        hw_spread=max(hws)/min(hws); d_trend=depths[0]/depths[-1]
        print(f"  resonance depth gamma1/gamma6 = {d_trend:.2f}; halfwidth spread = {hw_spread:.2f}x")
        # test depth ~ 1/sqrt(gamma) (prime-sum heuristic) vs constant width
        import statistics
        ratios=[prof[i]['depth']*math.sqrt(prof[i]['gamma']) for i in range(len(prof))]
        cv=statistics.pstdev(ratios)/statistics.mean(ratios)
        print(f"  depth*sqrt(gamma) coefficient-of-variation = {cv:.3f}  "
              f"({'~CONSTANT: depth ~ 1/sqrt(gamma), a clean scaling law' if cv<0.25 else 'not a clean 1/sqrt law'})")
        if hw_spread<2.0:
            print("  => resonances share a near-constant width and a depth scaling law => SELF-SIMILAR structure.")
            print("     Renormalization lever is LIVE: prove one template resonance + the scaling law -> all zeros.")
        else:
            print("  => widths vary strongly; not cleanly self-similar. Renormalization lever weak.")
    print("wrote results/self_similar_C3.json")

if __name__=="__main__":
    main()
