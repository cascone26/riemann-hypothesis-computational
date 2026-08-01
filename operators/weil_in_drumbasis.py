import numpy as np, os
import mpmath as mp

# LANE 1: the Weil / explicit-formula quadratic form in the EXPLICIT drum basis.
# Weil's criterion: Q(g) = sum_gamma |g^(gamma)|^2 - 2 sum_{p,k} (log p) p^{-k/2} g(k log p) g-real...
# In operator form with the zero-drum J (spectral measure = sum delta_gamma), for a test function h
# on the spectrum, sum_gamma h(gamma) = tr h(J). The Weil form's PRIME side is the explicit-formula
# dual. Concretely we test the a=1 RAZOR in the drum representation: build the finite Weil-form Gram
# on Gaussian bump test functions centered in [G0,G1], compare min-eig to the sine-Galerkin lambda_a.
# If the drum representation reproduces the razor (min-eig ~ 0), the RH-equivalent form is expressed
# in our explicit basis.

mp.mp.dps = 15

def archimedean_kernel(u):
    # the archimedean contribution to the explicit formula (digamma side), even in u
    # W_arch(u) = (1/2pi) Re psi(1/4 + i u/2) - (1/2) log pi  (per Weil/Suzuki normalization)
    return float(mp.re(mp.digamma(mp.mpf(0.25) + 1j * u / 2))) / (2 * np.pi) - 0.5 * np.log(np.pi) / (2 * np.pi)

def vonmangoldt(NM):
    L = np.zeros(NM + 1)
    for p in range(2, NM + 1):
        if L[p] == 0 and all(p % d for d in range(2, int(p ** 0.5) + 1)):
            pk = p
            while pk <= NM:
                L[pk] = np.log(p); pk *= p
    return L

# Weil quadratic form on test functions g_j (bumps): Q_{ij} = explicit-formula pairing.
# Use the standard smooth form: g_j(x)=exp(-(x-c_j)^2/(2 w^2)); its cosine transform hits zeros/primes.
# Q(g) = sum_gamma ghat(gamma)^2 - prime side - arch side. Build Gram over a set of c_j.
def ghat(g_c, g_w, t):
    # Fourier cosine transform of exp(-(x-c)^2/2w^2)+exp(-(x+c)^2/2w^2) at frequency t
    return 2 * g_w * np.sqrt(2 * np.pi) * np.exp(-0.5 * g_w ** 2 * t ** 2) * np.cos(c := g_c, ) * 0 + \
           2 * np.sqrt(2 * np.pi) * g_w * np.exp(-0.5 * g_w ** 2 * t ** 2) * np.cos(g_c * t)

DIR = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(DIR, "zeros10k.npy"))
gam = np.array([g for g in Z if g < 200])
NM = 400
L = vonmangoldt(NM)
prime_terms = [(n, L[n]) for n in range(2, NM + 1) if L[n] > 0]

cs = np.linspace(0.5, 4.0, 22)   # test-function centers in "x" (log-scale conjugate)
w = 0.35
K = len(cs)
Q = np.zeros((K, K))
for i in range(K):
    for j in range(K):
        # zero side: sum_gamma ghat_i(gamma) ghat_j(gamma)
        gi = 2 * np.sqrt(2 * np.pi) * w * np.exp(-0.5 * w ** 2 * gam ** 2) * np.cos(cs[i] * gam)
        gj = 2 * np.sqrt(2 * np.pi) * w * np.exp(-0.5 * w ** 2 * gam ** 2) * np.cos(cs[j] * gam)
        zero_side = np.sum(gi * gj)
        # prime side: sum_{n} L(n)/sqrt(n) * [g_i(log n) g_j(log n)] (even test fn value at x=log n)
        ps = 0.0
        for n, Ln in prime_terms:
            x = np.log(n)
            gival = np.exp(-0.5 * ((x - cs[i]) / w) ** 2) + np.exp(-0.5 * ((x + cs[i]) / w) ** 2)
            gjval = np.exp(-0.5 * ((x - cs[j]) / w) ** 2) + np.exp(-0.5 * ((x + cs[j]) / w) ** 2)
            ps += (Ln / np.sqrt(n)) * gival * gjval
        Q[i, j] = zero_side - ps
Q = 0.5 * (Q + Q.T)
ev = np.linalg.eigvalsh(Q)
print(f"Weil-form Gram (drum/explicit representation), {K} bump test fns, w={w}:")
print(f"  eigenvalues range [{ev[0]:.4e}, {ev[-1]:.4e}]")
print(f"  min eigenvalue = {ev[0]:.6e}  (a=1 razor: expect ~0^+ under RH; negative = truncation/prime-cutoff)")
print(f"  n_positive={np.sum(ev>1e-9*ev[-1])}, n_near_zero={np.sum(np.abs(ev)<1e-9*ev[-1])}")
print("bank: does the explicit-basis Weil form reproduce the razor (min-eig ~ 0)?")
