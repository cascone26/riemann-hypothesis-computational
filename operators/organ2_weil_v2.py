"""ORGAN 2 v2: get the Weil form RIGHT via the explicit-formula IDENTITY self-test.
Test function h even/real. Explicit formula (standard convention):
  sum_rho hhat(gamma_rho) = hhat(i/2)+hhat(-i/2) - sum_{n>=2} Lambda(n)/sqrt(n)[g(log n)+g(-log n)]
                            + (1/2pi) int hhat(r) [Re psi(1/4+ir/2) - log pi] dr
where g = inverse FT of hhat = h itself (g(x)=h(x)); hhat(r)=int h(x)e^{-irx}dx.
VALIDATE: LHS (sum over first N real zeros gamma) == RHS, for a chosen h. Only then build Q."""
import numpy as np, mpmath as mp
mp.mp.dps = 20
from sympy import primerange
z = np.load("zeros100k.npy")  # actual zeros (real, on line)
# h(x)=exp(-x^2/(2 t^2)) cos(a x) (even,real). g=h. hhat(r)= t sqrt(2pi)/2 [G(r-a)+G(r+a)], G=exp(-t^2 y^2/2)
t, a = 1.0, 3.0
def hhat(r): return t*np.sqrt(2*np.pi)/2*(np.exp(-t**2*(r-a)**2/2)+np.exp(-t**2*(r+a)**2/2))
def h(x): return np.exp(-x**2/(2*t**2))*np.cos(a*x)
# LHS: sum over zeros (+/- gamma) of hhat(gamma). Zeros come in +/-gamma pairs; hhat even so 2*sum_{gamma>0}
LHS = 2*np.sum(hhat(z))
# RHS terms:
# pole: hhat(i/2)+hhat(-i/2) = 2 Re hhat(i/2)  (hhat even)
def hhat_c(r): return t*mp.sqrt(2*mp.pi)/2*(mp.e**(-t**2*(r-a)**2/2)+mp.e**(-t**2*(r+a)**2/2))
pole = 2*float(mp.re(hhat_c(mp.mpf('0.5')*1j)))
# prime: sum_{p^m} log p / p^{m/2} [h(log p^m)+h(-log p^m)] = 2 sum log p/p^{m/2} h(log p^m) (h even)
prime = 0.0
for p in primerange(2,100000):
    lp=np.log(p); m=1
    while m*lp<20:
        prime += 2*lp/np.sqrt(p**m)*h(m*lp); m+=1
# arch: (1/2pi) int hhat(r)[Re psi(1/4+ir/2)-log pi] dr
rs=np.linspace(-50,50,8000)
Om=np.array([float(mp.re(mp.digamma(mp.mpf('0.25')+1j*r/2))-mp.log(mp.pi)) for r in rs])
arch = np.trapezoid(hhat(rs)*Om, rs)/(2*np.pi)
RHS = pole - prime + arch
print("EXPLICIT-FORMULA IDENTITY self-test (h=Gauss*cos, a=3):")
print(f"  LHS (sum over 2*100k zeros) = {LHS:.5f}")
print(f"  RHS = pole({pole:.4f}) - prime({prime:.4f}) + arch({arch:.4f}) = {RHS:.5f}")
print(f"  match? |LHS-RHS| = {abs(LHS-RHS):.5f}  (should be ~0 if formula correct; residual = zeros beyond 100k)")
print("DONE")
