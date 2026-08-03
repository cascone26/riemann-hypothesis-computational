import numpy as np, os

# Does the phase-only comb GENERALIZE to L(s,chi_3)? 1/L = sum mu(n)chi(n) n^{-s}, chi completely mult.
# Test G_L(u)=sum_gamma e^{-i arg L'(rho)} e^{iu gamma} vs mu(n)chi(n)/sqrt(n) on squarefree n.
DIR = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(DIR, "Lchi3_zeros.npy"))
aLp = np.load(os.path.join(DIR, "Lchi3_argLp.npy"))
print(f"L(s,chi_3): {len(Z)} zeros up to {Z[-1]:.1f}")

def factor(n):
    f = {}; d = 2; m = n
    while d * d <= m:
        while m % d == 0: f[d] = f.get(d, 0) + 1; m //= d
        d += 1
    if m > 1: f[m] = f.get(m, 0) + 1
    return f
def mu(n):
    f = factor(n)
    return 0 if any(e > 1 for e in f.values()) else (-1 if len(f) % 2 else 1)
chiv = {1: 1, 2: -1, 0: 0}
def chi(n):  # completely multiplicative mod 3
    r = 1
    for _ in range(1): pass
    return chiv[n % 3]

# mu_chi(n) = mu(n) * chi(n) (chi completely mult, so chi(n)=prod chi(p))
sf = [n for n in range(2, 1500) if mu(n) != 0 and chi(n) != 0]  # coprime to 3, squarefree
mus = np.array([mu(n) * chi(n) for n in sf], float)
absmus = np.abs(mus)
arr = np.array(sf, float)
E = np.exp(1j * np.outer(np.log(arr), Z))
w = np.exp(-1j * aLp)

def best(target):
    G = (E @ w) / np.sum(np.abs(w))
    return max(np.corrcoef(np.real(np.exp(1j * c) * G), target)[0, 1] for c in np.linspace(-np.pi, np.pi, 361))

cs_signed = best(mus / np.sqrt(arr))
cs_abs = best(absmus / np.sqrt(arr))
# control: random-phase weight
rng = np.random.default_rng(0)
ctrl = [best_r for best_r in [max(np.corrcoef(np.real(np.exp(1j*c)*((E@np.exp(1j*rng.uniform(-np.pi,np.pi,len(Z))))/len(Z))), mus/np.sqrt(arr))[0,1] for c in np.linspace(-np.pi,np.pi,121)) for _ in range(100)]]
ctrl = np.array(ctrl)
print(f"phase-only e^-i argL': corr vs SIGNED mu(n)chi(n)/sqrt(n) = {cs_signed:+.3f}")
print(f"                       corr vs UNSIGNED |mu chi|/sqrt(n)   = {cs_abs:+.3f}")
print(f"random-phase control: mean {ctrl.mean():+.3f} sd {ctrl.std():.3f} max {ctrl.max():+.3f}")
print(f"z = {(cs_signed - ctrl.mean())/ctrl.std():+.2f}")
print("\nif signed >> control: the phase-only Mobius comb GENERALIZES to L-functions (chi-twisted Mobius)")
print("=> it is a general Selberg-class phenomenon, not special to zeta. Strengthens novelty massively.")
