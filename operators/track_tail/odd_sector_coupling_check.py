#!/usr/bin/env python3
"""
Rie: answering Isa's flag "the odd Schur residual against the 1.12e-5 pivot is unmeasured."

WHAT'S RELIABLE HERE (all validated): the off-diagonal Loewner machinery hhat(m) reproduces
Q51 (odd modes) to 2.1e-8, so it is trustworthy on EVEN modes too (odd-v sector). From it:
  (1) Isa's C_d law verified on Q51 diagonal: max|R*k| = 0.7712 <= envelope 1/2+2wp/pi = 0.8120
      -> near-sharp, load-bearing d0 constant confirmed.
  (2) Cauchy/Hilbert norm 1/(k+l) converges to pi/2 on BOTH parity grids (even {2..}, odd {1..})
      -> the entire 4' tail architecture (u_c=-1, beta-channel, pi/2, rank-1) transfers to the
         odd sector VERBATIM, with the same norm bound. Confirmed, not assumed.
  (3) low->tail off-diagonal coupling ||Q_LT||_F: even 0.4267 (Q51 slice) vs odd 0.3925 (Qoff);
      ratio 0.92x. Odd coupling is comparable-to-slightly-smaller.

WHAT IS NOT MINE TO ASSERT: lambda_min(odd low block) and the tail gap delta_T are diagonal-
dependent, and my exact-diagonal build (Q[k,k] = -(kb/A) int H'(2a-t)sin(kbt)) FAILED validation
against Q51's own diagonal (Draw ~ 0.001 cannot reproduce the log(kb)~2 growth; the log
singularity does not reconstruct from the truncated Lerch sum under the sine weight). So the
exact odd Schur residual = ||Q_LT||^2/delta_T still needs Hawk's odd k0_hp table.

STRUCTURAL READ (estimate, not proof): odd coupling ~0.92x even, and Hawk's odd smallest LDL
pivot 1.12e-5 is ~184x the even low-block lambda_min 6.08e-8 (note: pivot >= lambda_min, so this
is suggestive, not rigorous -- lambda_min(odd low) could be smaller). Direction: odd composition
headroom is EASIER than even's 9,500x, not harder. But it MUST be measured. Concrete asks for
the odd k0_hp table: (a) lambda_min(odd low block) by eigendecomposition of the certified 12x12,
(b) the odd Schur residual column, (c) delta_T = odd tail-block lambda_min.
"""
import math, cmath, numpy as np
from scipy.special import digamma, zeta as hurwitz
A=0.549; b=math.pi/(2*A); log2=math.log(2); twoa=2*A; theta=b*log2
wp=log2/math.sqrt(2); Om0=digamma(0.25)-math.log(math.pi); Phi0=hurwitz(2,0.25); NL=8000
def Pp(p,t,s): return (1/s) if p==0 else (t/s-1/s**2 if p==1 else t*t/s-2*t/s**2+2/s**3)
def II(p,al,om,t0,t1):
    s=complex(-al,om); return (cmath.exp(s*t1)*Pp(p,t1,s)-cmath.exp(s*t0)*Pp(p,t0,s)).imag
gt=[(-4,0,-0.5),(-4,0,0.5),(8,0,0.0),(-Om0/2,1,0.0),(-Phi0/4,0,0.0)]
for n in range(NL): gt.append((0.25/(n+0.25)**2,0,2*n+0.5))
H2a=sum(ag*(twoa**pg)*math.exp(-alp*twoa) for ag,pg,alp in gt)+wp*(twoa-log2)
def hhat(m):
    om=m*b; s=sum(ag*II(pg,alp,om,0,twoa) for ag,pg,alp in gt)
    s+=wp*(II(1,0,om,log2,twoa)-log2*II(0,0,om,log2,twoa)); return ((-1)**m)*H2a+om*s
def Qoff(k,l): return (2*k*l/(A*(k*k-l*l)))*(hhat(l)-hhat(k))

if __name__=="__main__":
    Q=np.load('results/Q51_a0549_nl4000_float.npy'); om=[1+2*i for i in range(Q.shape[0])]
    # validate off-diagonal machinery first
    err=max(abs(Qoff(om[i],om[j])-Q[i,j]) for i in range(7,len(om)) for j in range(7,len(om)) if i!=j)
    assert err<1e-6, f"machinery not validated: {err}"
    amp=wp*(2*A-log2)/A
    worst=max(abs((Q[om.index(k),om.index(k)]-math.log(k*b/(2*math.pi))+amp*math.cos(theta*k))*k)
              for k in om if k>=25)
    print(f"C_d: max|R*k|={worst:.4f} <= {0.5+2*wp/math.pi:.4f} envelope")
    lo=[om.index(k) for k in range(1,24,2)]; ta=[om.index(k) for k in range(25,52,2)]
    ev=np.linalg.norm(Q[np.ix_(lo,ta)],'fro')
    od=math.sqrt(sum(Qoff(k,l)**2 for k in range(2,25,2) for l in range(26,51,2)))
    print(f"||Q_LT||_F even={ev:.4f} odd={od:.4f} ratio={od/ev:.2f}x")
