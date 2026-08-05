"""UNCONDITIONAL Weil quadratic form W(phi) = arch + pole - prime, RH <=> W>=0 for all phi.
NO zeros used - all three terms explicit. Map the positivity boundary; find where prime nearly
cancels arch+pole (the extremal, 'barely true' region). Test function family: Gaussian-modulated
cosines phi_hat(r) concentrated at frequency mu, width sigma."""
import numpy as np, mpmath as mp
mp.mp.dps = 20
from sympy import primerange
# Weil explicit formula for h(r)=|phihat(r)|^2 (h>=0), with test phi. Use the standard form:
# sum_rho h(gamma) = (1/2pi) int_{-inf}^{inf} h(r) Omega(r) dr + h(i/2)+h(-i/2) - sum_n Lambda(n)/sqrt(n) ghat(log n)*2
# Omega(r) = Re[psi(1/4 + i r/2)] - log(pi)   (archimedean density)
# ghat = fourier transform; for h(r)=e^{-(r-mu)^2/(2s^2)} type... use real even test functions.
# Simplify: take g(x) (even, real, in x=log-variable), h(r)=ghat(r)=int g(x)e^{-irx}dx.
# W = arch + pole - prime where:
#   arch = (1/2pi) int h(r) Omega(r) dr
#   pole = h(i/2)+h(-i/2)  [from the s=1 pole; = 2*analytic-continuation of h to r=i/2]
#   prime = 2 sum_n Lambda(n)/sqrt(n) g(log n)
# Use g(x) = exp(-x^2/(2 t^2)) cos(a x)  (real even), h(r) = ghat(r).
def Omega(r):
    return float(mp.re(mp.digamma(mp.mpf('0.25') + 1j*r/2)) - mp.log(mp.pi))
def run(tw, a):
    # g(x)=exp(-x^2/(2 tw^2)) cos(a x); ghat(r)= tw*sqrt(2pi)/2 [e^{-tw^2(r-a)^2/2}+e^{-tw^2(r+a)^2/2}]
    def h(r): return tw*np.sqrt(2*np.pi)/2*(np.exp(-tw**2*(r-a)**2/2)+np.exp(-tw**2*(r+a)**2/2))
    rs=np.linspace(-60/max(tw,0.1),60/max(tw,0.1),8000)
    arch=np.trapezoid(h(rs)*np.array([Omega(r) for r in rs]),rs)/(2*np.pi)
    # pole: h analytically continued to r=i/2: replace (r-a)->(i/2-a) etc.
    pole=2*float(mp.re(tw*mp.sqrt(2*mp.pi)/2*(mp.e**(-tw**2*(1j/2-a)**2/2)+mp.e**(-tw**2*(1j/2+a)**2/2))))
    prime=0.0
    for p in primerange(2,int(np.exp(60/max(tw,0.1)))+2 if tw>0.3 else 100000):
        lp=np.log(p); m=1
        while m*lp<12*tw+8:
            g_val=np.exp(-(m*lp)**2/(2*tw**2))*np.cos(a*m*lp)
            prime+=2*lp/np.sqrt(p**m)*g_val
            m+=1
        if p>50000: break
    return arch, pole, prime, arch+pole-prime
print(" width t  freq a   arch     pole     prime     W=arch+pole-prime  (>=0 unconditionally?)")
for tw in (2.0, 1.0, 0.5):
    for a in (0.0, 2.0, 5.0, 10.0):
        arch,pole,prime,W=run(tw,a)
        print(f"  {tw:.1f}    {a:4.1f}   {arch:+.3f}  {pole:+.3f}  {prime:+.3f}   {W:+.4f}  {'POS' if W>=0 else 'NEG!'}", flush=True)
print("DONE")
