"""ORGAN 2 v1: Weil quadratic form as a matrix on a test-function basis -> SDP/SOS lane.
Polarized bilinear form Q_{jk}=W(phi_j,phi_k)=arch+pole-prime. RH <=> Q >= 0 (PSD).
Basis: phi_k(x)=exp(-x^2/(2 t^2)) cos(a_k x). Build Q, check PSD (RH-consistency check),
then this becomes the object we certify (exact rational Gram) + run the hierarchy on."""
import numpy as np, mpmath as mp
mp.mp.dps = 25
from sympy import primerange
def Omega(r): return float(mp.re(mp.digamma(mp.mpf('0.25')+1j*r/2))-mp.log(mp.pi))
# For phi with phihat(r), the Weil bilinear form entries (h=phihat_j * conj(phihat_k), even real):
# arch = (1/2pi) int phihat_j(r) phihat_k(r) Omega(r) dr
# pole = phihat_j(i/2)phihat_k(i/2)+... (analytic continuation) [approx via the h(i/2) terms]
# prime = 2 sum_n Lambda(n)/sqrt(n) g_j*g_k-ish -> use g_jk(x)=phi_j conv... simplest: work in x-space:
#   prime term for the PRODUCT: 2 sum_n Lambda(n)/sqrt(n) (phi_j * phi_k~)(log n).
# Use phi_k(x)=exp(-x^2/(2t^2))cos(a_k x); phihat_k(r)= (t sqrt(2pi)/2)[G(r-a_k)+G(r+a_k)], G(y)=exp(-t^2 y^2/2)
t = 1.0
avals = [0.0, 4.0, 8.0, 12.0]
n = len(avals)
def phihat(k, r):
    a=avals[k]
    return t*np.sqrt(2*np.pi)/2*(np.exp(-t**2*(r-a)**2/2)+np.exp(-t**2*(r+a)**2/2))
def phix(k, x):  # phi_k in x-space
    return np.exp(-x**2/(2*t**2))*np.cos(avals[k]*x)
rs = np.linspace(-40,40,6000)
Om = np.array([Omega(r) for r in rs])
primes_pp = []
for p in primerange(2, 100000):
    lp=np.log(p); m=1
    while m*lp < 30:
        primes_pp.append((m*lp, lp/np.sqrt(p**m))); m+=1
Q = np.zeros((n,n))
for j in range(n):
    for k in range(n):
        arch = np.trapezoid(phihat(j,rs)*phihat(k,rs)*Om, rs)/(2*np.pi)
        # pole: phihat(i/2) via analytic continuation (r->i/2): G(i/2-a)=exp(-t^2(i/2-a)^2/2)
        def phihat_c(kk):
            a=avals[kk]
            return t*np.sqrt(2*np.pi)/2*(mp.e**(-t**2*(mp.mpf('0.5')*1j-a)**2/2)+mp.e**(-t**2*(mp.mpf('0.5')*1j+a)**2/2))
        pole = float(mp.re(phihat_c(j)*mp.conj(phihat_c(k))))
        # prime: 2 sum (phi_j * phi_k~)(log p^m) ~ 2 sum Lambda/sqrt(n) phi_j(L)phi_k(L) (diagonal approx)
        prime = 2*sum(w*phix(j,L)*phix(k,L) for L,w in primes_pp)
        Q[j,k] = arch + pole - prime
Q = (Q+Q.T)/2
ev = np.linalg.eigvalsh(Q)
print("Weil form matrix Q (arch+pole-prime) on 4 Gaussian test functions:")
print("  eigenvalues:", np.array2string(ev, precision=4))
print("  min-eig = %.5f  => %s (RH predicts PSD, min-eig>=0)" % (ev[0], "PSD (RH-consistent)" if ev[0]>-1e-9 else "NEG - check formula!"))
print("  Q symmetric:", np.allclose(Q,Q.T))
print("DONE")
