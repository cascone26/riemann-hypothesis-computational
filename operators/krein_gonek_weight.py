import numpy as np, os
import mpmath as mp

# Sharpen §517: compare three zero-sum weights on the mu(n)/sqrt(n) target:
#   W_phase:  e^{-2i arg zeta'(rho)}          (unit modulus — the string's field)
#   W_gonek:  1/zeta'(rho)                    (Perron/Gonek exact weight)
#   W_half:   e^{-i arg zeta'(rho)}           (single twist control)
# plus plain (no weight) as anti-control. Metric: best-rotation corr with mu(n)/sqrt(n).

mp.mp.dps = 16
DIR = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(DIR, "zeros500.npy"))
ph = np.load(os.path.join(DIR, "argzp500.npy"))
fn = os.path.join(DIR, "abszp500.npy")
if os.path.exists(fn):
    mag = np.load(fn)
else:
    mag = np.empty(len(Z))
    for i, g in enumerate(Z):
        mag[i] = float(abs(mp.zeta(mp.mpc(0.5, float(g)), derivative=1)))
        if i % 100 == 0:
            print(f"  |zeta'| {i}/{len(Z)}", flush=True)
    np.save(fn, mag)

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

sf = [n for n in range(2, 99) if mu(n) != 0]
target = np.array([mu(n) / np.sqrt(n) for n in sf])
lgn = np.log(np.array(sf, dtype=float))
E = np.exp(1j * np.outer(lgn, Z))

def score(w, name):
    Gn = (E * w[None, :]).sum(axis=1) / np.sum(np.abs(w))
    best = max((float(np.corrcoef(np.real(np.exp(1j * c) * Gn), target)[0, 1]), float(c))
               for c in np.linspace(-np.pi, np.pi, 721))
    print(f"  {name:>28}: best r = {best[0]:+.4f} (c={best[1]:+.2f})")
    return best[0]

print(f"{len(sf)} squarefree targets, 500 zeros:")
score(np.ones(len(Z)), "plain (no weight)")
score(np.exp(-1j * ph), "e^{-i arg zeta'}")
score(np.exp(-2j * ph), "e^{-2i arg zeta'} (string)")
score(np.exp(-2j * ph) / mag, "1/|z'| * e^{-2i arg} ")
score(np.exp(-1j * ph) / mag, "1/zeta'(rho) (Gonek)")
score(1.0 / mag * np.ones(len(Z)) * np.exp(0j), "1/|zeta'| (mag only)")
