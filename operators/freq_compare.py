#!/usr/bin/env python3
"""§898: FREQUENCY-DOMAIN arch-vs-prime (Guinand–Weil frame). From §897, Ω(r)=C−Σ_n a_n/(a_n²+r²/4), a_n=n+¼.
Cosine transform of a_n/(a_n²+r²/4) is π e^{−2a_n|u|}, so the archimedean weight (u≠0 part) has CLOSED FORM
  Ŵ(u) := cosine-transform of Ω = −π Σ_{n≥0} e^{−2(n+¼)|u|} = −π e^{−|u|/2}/(1−e^{−2|u|}).
VERIFY this closed form numerically. Then frame the comparison: the arch form in frequency is a Toeplitz/convolution
operator with kernel Ŵ(u−u')+Ŵ(u+u') (+ a positive diagonal from C), the prime form is point masses Λ(n)/√n at
u=log n. This is the natural place the domination arch⪰prime lives. Report the closed form + the prime masses in the
band |u|<4W and their relative scale."""
import numpy as np
from scipy.special import digamma
np.set_printoptions(precision=5,suppress=True,linewidth=140)

def Wclosed(u):
    u=np.abs(u); out=np.empty_like(u)
    small=u<1e-8
    out[~small]=-np.pi*np.exp(-u[~small]/2)/(1-np.exp(-2*u[~small]))
    out[small]=-np.inf  # diverges at 0 (the C·δ sits there separately)
    return out

# numeric cosine transform of Omega(r) over a tapered range, minus the divergent mean (compare shape to closed form)
def Wnumeric(u, R=4000.0, dr=0.01):
    r=np.arange(dr,R,dr); Om=np.real(digamma(0.25+0.5j*r))-np.log(np.pi)
    # subtract large-r asymptote ln(r/2pi) so the transform converges (that piece is the u->0 singularity/const)
    asymp=np.log(r/(2*np.pi))
    g=Om-asymp
    taper=0.5*(1+np.cos(np.pi*r/R))  # cosine taper to R
    vals=[]
    for uu in u:
        vals.append(2*np.sum(g*np.cos(uu*r)*taper)*dr)
    return np.array(vals)

print("="*100)
print("§898 archimedean weight Ŵ(u): closed form −π e^{−u/2}/(1−e^{−2u}) vs numeric cosine transform of (Ω − ln(r/2π))")
print("="*100)
us=np.array([0.3,0.5,0.6931,1.0,1.0986,1.5,2.0,2.5,3.0,3.5])
Wc=Wclosed(us); Wn=Wnumeric(us)
print(f"{'u':>8} {'closed Ŵ(u)':>13} {'numeric':>11}  {'note':>16}")
for i,uu in enumerate(us):
    note=""
    for n in range(2,40):
        if abs(uu-np.log(n))<0.01: note=f"= log {n}"
    print(f"{uu:>8.4f} {Wc[i]:>13.4f} {Wn[i]:>11.4f}  {note:>16}")
print("-"*100)
# prime masses in band |u|<4 (4W=4): Lambda(n)/sqrt(n) at u=log n
print("prime masses Λ(n)/√n at u=log n (|u|<4, i.e. n<e^4≈54.6):")
def vm(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
LAM=vm(60)
for n in range(2,55):
    if LAM[n]>0:
        u=np.log(n); wc=float(Wclosed(np.array([u]))[0]); mass=LAM[n]/np.sqrt(n)
        print(f"  n={n:2d} u=log n={u:.4f}  Λ/√n={mass:.5f}   |Ŵ(u)|(arch weight)={-wc:.4f}   ratio arch/prime={-wc/mass:.2f}")
print("="*100)
print("The archimedean weight |Ŵ(u)| vs the prime mass Λ(n)/√n at each u=log n: where |Ŵ|>Λ/√n the arch dominates")
print("pointwise in frequency. This is the Guinand–Weil-frame version of arch⪰prime; the ratio's behavior is the lever.")
print("="*100)
