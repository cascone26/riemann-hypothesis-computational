#!/usr/bin/env python3
"""§199.6: CONTINUUM via a-INTERVAL Arb. Certify the even 4×4 PD over a PANEL a∈[a_lo,a_hi] in one shot (A is an Arb
interval; all closed-form entries propagate intervals; term-by-term Lerch). Only p=2 contributes on the whole window
(2a∈(0.84,1.099)⊂(log2,log3)). Test one panel; then sweep (0.42,0.5493)."""
import sys
from flint import arb, acb, ctx
ctx.prec=1600
NLERCH=8000
def certify_panel(alo, ahi):
    A=arb(str(alo)).union(arb(str(ahi)))     # interval [alo,ahi]
    B=arb.pi()/(2*A); LOG2=arb(2).log(); TWOA=2*A
    OM0=(arb(1)/4).digamma()-arb.pi().log(); PHI0=acb(2).zeta(acb(arb(1)/4)).real; WP=arb(2).log()/arb(2).sqrt()
    def Pp(p,t,s):
        if p==0: return 1/s
        if p==1: return t/s-1/s**2
        return t*t/s-2*t/s**2+2/s**3
    def II(p,alpha,omega,sc,t0,t1):
        s=acb(-alpha,omega); F=lambda tt:(s*acb(tt)).exp()*Pp(p,acb(tt),s); v=F(t1)-F(t0)
        return v.real if sc=='cos' else v.imag
    def dterms(k,l):
        if k==l: return [(2*A*arb(k)**2*B**2,0,arb(k)*B,'cos'),(-arb(k)**2*B**2,1,arb(k)*B,'cos'),(-arb(k)*B,0,arb(k)*B,'sin')]
        al=arb(k)*arb(l)*B*(arb(1)/(k-l)-arb(1)/(k+l)); ak=-arb(k)*arb(l)*B*(arb(1)/(k-l)+arb(1)/(k+l))
        return [(al,0,arb(l)*B,'sin'),(ak,0,arb(k)*B,'sin')]
    def Dmax(k,l):
        tot=arb(0)
        for (ad,pd,om,sc) in dterms(k,l): tot+=abs(ad)*(TWOA if pd==1 else arb(1))
        return tot
    def Q(k,l):
        tot=arb(0)
        gt=[(arb(-4),0,arb(-1)/2),(arb(-4),0,arb(1)/2),(arb(8),0,arb(0)),(-OM0/2,1,arb(0)),(-PHI0/4,0,arb(0))]
        for n in range(NLERCH): gt.append((arb(1)/4/(arb(n)+arb(1)/4)**2,0,2*arb(n)+arb(1)/2))
        for (ag,pg,alp) in gt:
            for (ad,pd,om,sc) in dterms(k,l): tot+=ag*ad*II(pg+pd,alp,om,sc,arb(0),TWOA)
        tb=Dmax(k,l)/8*(arb(1)/(2*(arb(NLERCH)-1)**2)); tot+=tb*arb(0,1)
        for (ad,pd,om,sc) in dterms(k,l):
            tot+=WP*ad*II(pd+1,arb(0),om,sc,LOG2,TWOA)+(-WP*LOG2)*ad*II(pd,arb(0),om,sc,LOG2,TWOA)
        return tot
    modes=[5,7,3,1]; M=[[Q(k,l) for l in modes] for k in modes]
    d=[None]*4; L=[[arb(0)]*4 for _ in range(4)]
    for i in range(4):
        s=M[i][i]
        for kk in range(i): s=s-L[i][kk]**2*d[kk]
        d[i]=s
        if not (s>0): return False, i, s.str(4)
        for j in range(i+1,4):
            t=M[j][i]
            for kk in range(i): t=t-L[j][kk]*L[i][kk]*d[kk]
            L[j][i]=t/d[i]
    return True, None, d[3].str(4)
if __name__=="__main__":
    alo,ahi=(float(sys.argv[1]),float(sys.argv[2])) if len(sys.argv)>2 else (0.45,0.46)
    ok,fp,info=certify_panel(alo,ahi)
    print(f"§199.6 panel a∈[{alo},{ahi}] (width {ahi-alo}): {'CERTIFIED PD' if ok else f'FAIL at pivot {fp}'}  (last pivot {info})")
