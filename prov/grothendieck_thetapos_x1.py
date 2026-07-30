import numpy as np
from scipy.integrate import solve_ivp
# TARGET (theta_pos lane): the x=1 (finite, "local/archimedean") endpoint boundary form.
# L u = -((1-x^2)u')' + (2pi)^2 x^2 u = mu u,  mu=-E^2,  p(x)=1-x^2, on (1,inf).
# Energy quadratic form:  <Lu,u> = INT[ (1-x^2)|u'|^2 + (2pi)^2 x^2 |u|^2 ] - [ (1-x^2) u' u ]_1^X.
# The x=1 boundary term  B1 = lim_{x->1+} (1-x^2) u'(x) u(x)  is LOCAL (X-independent).
# Extension at x=1: u = cos(b) u_r + sin(b) u_l,  u_r~1 (regular), u_l~log(x-1) (singular). b=0 = Friedrichs = build.
# PRE-REGISTERED (ink, before run):
#   (P1) B1 diverges like C(b)*log(1/eps) as eps->0 for b!=0, coeff C(b) ~ sin^2(b) (=> energy UNbounded below).
#   (P2) At b=0 (Friedrichs) the divergent coeff C(0)=0 within tol -> energy bounded below (positive extension).
#   => theta_pos at the x=1 factor is the ISOLATED point b=0, and it CONTAINS the build base point B(E)=0.
FP=(2*np.pi)**2
def B1_vs_eps(E,beta,X=40.0,epss=(1e-3,3e-4,1e-4,3e-5,1e-5)):
    mu=-E*E
    def rhs(x,y): return [y[1],(2*x*y[1]+(FP*x*x-mu)*y[0])/(1-x*x)]
    out=[]
    for eps in epss:
        ICr=np.array([1.0,(mu-FP)/2.0])            # regular:  u_r(1)=1
        ICl=np.array([np.log(eps),1.0/eps])        # singular: u_l ~ log(x-1)
        IC=np.cos(beta)*ICr+np.sin(beta)*ICl
        x0=1+eps
        B1=(1-x0*x0)*IC[1]*IC[0]                    # (1-x^2) u' u  AT the endpoint x=1+eps
        out.append((eps,B1))
    return out
gam=14.134725
print("theta_pos / x=1 boundary form  B1 = (1-x^2)u'u  vs eps, at E=gamma_1=%.5f\n"%gam)
print(" beta      B1(eps=1e-3)  1e-4    1e-5    | log-divergence coeff C(beta)=dB1/dlog(1/eps)")
for beta in [0.0,0.15,0.3,0.6,1.0,1.5]:
    d=B1_vs_eps(gam,beta)
    eps=np.array([e for e,_ in d]); B=np.array([b for _,b in d])
    # fit B1 ~ C*log(1/eps) + const
    C=np.polyfit(np.log(1/eps),B,1)[0]
    print(f" {beta:5.3f}   {B[0]:11.3f} {B[2]:8.2f} {B[4]:8.1f}  |  C={C:10.3f}   sin^2(b)={np.sin(beta)**2:.4f}")
print("\nRead: C(beta)=0 => energy boundary form finite (positive/lower-semibounded extension).")
print("      C(beta)<0 growing => energy -> -inf as eps->0 => NOT a positive extension.")
