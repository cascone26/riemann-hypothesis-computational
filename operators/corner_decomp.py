#!/usr/bin/env python3
"""§321: decompose the corner law. min gap is ALWAYS pair #0 (deterministic) — so min/mean = gap0/mean is a
DETERMINISTIC asymptotics question, not an RMT statistic (correcting §319's GUE reading if this holds).
Measure separately: t1, gap0=t2-t1, t_max, mean=(t_max-t1)/(d-1); local exponents of each; and test the saddle
prediction t_max ~ 2 d^{3/2}/u*(d) (from moment ratio gamma_{k}/gamma_{k+1} ~ (2k)^2/u*^2).
Algebra check: exp(min/mean) = exp(gap0) - exp(mean)."""
import mpmath as mp, pickle, math
mp.mp.dps=50
gam=[mp.mpf(g) for g in pickle.load(open("operators/_gam130.pkl","rb"))]
def usaddle(m):
    return mp.findroot(lambda u: 2*m/u+mp.mpf(9)/2-2*mp.pi*mp.e**(2*u), mp.mpf(0.5)+mp.log(m+1)/2)
rows=[]
print(" d :    t1       gap0      t_max     mean    min/mean   t_max_saddle_pred  ratio",flush=True)
for d in [16,25,36,50,70]:
    co=[mp.binomial(d,k)*gam[k] for k in range(d,-1,-1)]
    scale=mp.sqrt(abs(co[0])*abs(co[-1])); co=[c/scale for c in co]
    rts=mp.polyroots(co,maxsteps=500,extraprec=260)
    ts=sorted(mp.sqrt(abs(mp.re(r))) for r in rts)
    gap0=ts[1]-ts[0]; tmax=ts[-1]; mean=(ts[-1]-ts[0])/(d-1)
    us=usaddle(d)
    pred=2*d**mp.mpf(1.5)/us
    rows.append((d,ts[0],gap0,tmax,mean,gap0/mean))
    print(f"{d:3d}: {mp.nstr(ts[0],4):>7} {mp.nstr(gap0,4):>8} {mp.nstr(tmax,5):>9} {mp.nstr(mean,4):>8} {mp.nstr(gap0/mean,4):>8}   {mp.nstr(pred,5):>10}  {mp.nstr(tmax/pred,4)}",flush=True)
print("",flush=True)
print("local exponents:  t1     gap0    t_max   mean    min/mean  [check: gap0-mean]",flush=True)
for i in range(1,len(rows)):
    d0=rows[i-1]; d1=rows[i]
    l=lambda j: math.log(float(d1[j]/d0[j]))/math.log(d1[0]/d0[0])
    print(f"  d={d0[0]}->{d1[0]}: {l(1):+.3f} {l(2):+.3f} {l(3):+.3f} {l(4):+.3f}  {l(5):+.3f}   [{l(2)-l(4):+.3f}]",flush=True)
print("",flush=True)
print("interpretation: if t_max tracks the saddle law (ratio ~const) and gap0 has its own clean exponent,",flush=True)
print("the corner law is DETERMINISTIC moment-asymptotics (gap0 exp − t_max exp + 1), not an RMT statistic.",flush=True)
print("DONE",flush=True)
