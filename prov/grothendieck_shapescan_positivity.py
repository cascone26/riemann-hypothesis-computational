import numpy as np
from scipy.integrate import solve_ivp
# GUARD for the shape scan (theta_pos lane): does the POSITIVITY GEOMETRY survive varying the
# potential exponent b?  Generalized symbol:  L u = -((1-x^2)u')' + FP*x^(2b) u = mu u,  b=1 = build.
# Two things must hold for theta_pos = {beta=0} x (inf-circle) to stay well-posed across the scan:
#   (I)  x=1 endpoint: indicial structure {1,log} => C(beta)=2 sin^2 b  must be b-INDEPENDENT
#        (it's governed by p(x)=1-x^2's simple zero + q(1) finite, NOT by the exponent b).
#   (II) inf endpoint stays LIMIT-CIRCLE (both solns L^2) so the extension theta-circle survives.
#        WKB predicts u ~ x^-((b+1)/2), |u|^2 ~ x^-(b+1) => L^2 iff b>0.  Measure the decay exponent.
FP=(2*np.pi)**2
def C_beta_x1(E,beta,b,epss=(1e-3,1e-4,1e-5)):
    mu=-E*E
    def rhs(x,y): return [y[1],(2*x*y[1]+(FP*x**(2*b)-mu)*y[0])/(1-x*x)]
    B=[]
    for eps in epss:
        ICr=np.array([1.0,(mu-FP)/2.0]); ICl=np.array([np.log(eps),1.0/eps])
        IC=np.cos(beta)*ICr+np.sin(beta)*ICl; x0=1+eps
        B.append((1-x0*x0)*IC[1]*IC[0])
    eps=np.array(epss); return np.polyfit(np.log(1/eps),np.array(B),1)[0]
def tail_decay(E,b,X=60.0):
    mu=-E*E
    def rhs(x,y): return [y[1],(2*x*y[1]+(FP*x**(2*b)-mu)*y[0])/(1-x*x)]
    xs=np.linspace(20,X,4000)
    s=solve_ivp(rhs,[1+1e-4,X],[1.0,(mu-FP)/2.0],rtol=1e-10,atol=1e-12,max_step=0.02,t_eval=xs)
    u=s.y[0]
    # envelope: local maxima of |u|, fit log|env| vs log x
    a=np.abs(u); pk=(a[1:-1]>a[:-2])&(a[1:-1]>a[2:]); xp=xs[1:-1][pk]; ap=a[1:-1][pk]
    if len(xp)<4: return np.nan
    return np.polyfit(np.log(xp),np.log(ap),1)[0]   # decay exponent d: |u|~x^d
g=14.134725
print("(I)  x=1 divergence coeff C(beta) vs exponent b  [pre-reg: C=2sin^2(beta), b-INDEPENDENT]")
print("      b      C(b=0)   C(0.6)   C(1.0) | 2sin^2 ref: 0.0 / 0.638 / 1.416")
for b in [0.6,1.0,1.4,2.0]:
    c0=C_beta_x1(g,0.0,b); c6=C_beta_x1(g,0.6,b); c1=C_beta_x1(g,1.0,b)
    print(f"    {b:4.2f}   {c0:7.3f}  {c6:7.3f}  {c1:7.3f}")
print("\n(II) inf-endpoint decay exponent d (|u|~x^d)  [pre-reg: d=-(b+1)/2, L^2 iff b>0 => LC preserved]")
print("      b     measured d    predicted -(b+1)/2   L^2 at inf?")
for b in [0.4,0.7,1.0,1.4,2.0]:
    d=tail_decay(g,b); pred=-(b+1)/2
    L2 = "YES (LC)" if (2*d < -1) else "NO (LP)"
    print(f"    {b:4.2f}   {d:8.3f}      {pred:8.3f}          {L2}")
