"""
THE decisive step-(b) test: raw low-even-prolate span {ps_0,ps_2} FAILS to capture the Weil minimizer θ_x
(residual GROWS with window a: 0.009→0.21).  Connes' claim is k_λ = E(prolate): the ζ-multiplication map
E injects the prime structure.  In the additive log-variable x=log u, E(f)(x) = e^{x/2} Σ_{n≥1} f(x+log n)
(a comb of log-n shifts).  TEST: does span{E(ps_0), E(ps_2)} (restricted to [-a,a]) capture θ_x with a
residual that SHRINKS as a grows — the opposite of the raw-prolate failure?  If yes: Connes' k_λ mechanism
is confirmed on our own numerics.
"""
import numpy as np
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def weil_min(a, M):
    Nmax = int(np.floor(np.exp(2 * a))) + 1
    Lam = vonmangoldt_upto(min(Nmax, 30_000_000))
    x = np.linspace(-a, a, M + 1); hstep = x[1] - x[0]
    diffs = np.arange(-M, M + 1) * hstep
    gv = np.array([screw_g(d, Lam, Nmax) for d in diffs])
    idx = np.arange(M + 1); G = gv[(idx[:, None] - idx[None, :]) + M]
    w = np.sqrt(np.full(M + 1, hstep)); w[0] = w[-1] = np.sqrt(0.5 * hstep)
    Gs = w[:, None] * G * w[None, :]; Gs = (Gs + Gs.T) / 2
    ev, V = np.linalg.eigh(Gs); e = V[:, 0] / w
    return x, hstep, e / np.sqrt(np.sum(e * e) * hstep), ev[0]

def even_prolates(a, x, hstep, c_bw, k=2):
    D = x[:, None] - x[None, :]
    K = np.where(np.abs(D) < 1e-12, c_bw / np.pi, np.sin(c_bw * D) / (np.pi * D)) * hstep
    ev, V = np.linalg.eigh(K)
    evens = []
    for j in range(len(ev) - 1, -1, -1):
        v = V[:, j]
        if np.max(np.abs(v - v[::-1])) < 1e-6 * np.max(np.abs(v)):
            evens.append(v / np.sqrt(np.sum(v * v) * hstep))
            if len(evens) >= k:
                break
    return evens

def E_map(f, x, a, hstep, Nmax=None):
    """E(f)(x) = e^{x/2} sum_{n>=1} f(x+log n), f=0 outside [-a,a] (compact support)."""
    if Nmax is None:
        Nmax = int(np.floor(np.exp(2 * a))) + 2
    out = np.zeros_like(x)
    for n in range(1, Nmax + 1):
        xn = x + np.log(n)
        fv = np.interp(xn, x, f, left=0.0, right=0.0)  # zero outside window
        out += fv
        if np.log(n) > 2 * a:  # shift pushes support entirely out of [-a,a]
            break
    return np.exp(x / 2) * out

def norm(f, hstep):
    nn = np.sqrt(np.sum(f * f) * hstep)
    return f / nn if nn > 0 else f

def residual_in_span(theta, basis, hstep):
    """1 - ||P_span theta||^2 for orthonormalized real basis; theta normalized."""
    B = []
    for b in basis:
        v = b.copy()
        for u in B:
            v = v - (np.sum(v * u) * hstep) * u
        nn = np.sqrt(np.sum(v * v) * hstep)
        if nn > 1e-9:
            B.append(v / nn)
    proj2 = sum((np.sum(theta * u) * hstep) ** 2 for u in B)
    return 1.0 - proj2

print("Weil θ_x residual: raw {ps0,ps2}  vs  E-mapped {E ps0, E ps2}  (Connes k_λ), vs window a")
print(f"{'a':>5} {'ε(min)':>9} {'raw resid':>11} {'E-map resid':>12} {'better?':>8}")
for a in [1.0, 1.5, 2.0, 2.5, 3.0]:
    M = min(1400, max(700, int(900 / a)))
    x, hstep, theta, eps = weil_min(a, M)
    ps = even_prolates(a, x, hstep, c_bw=max(1.0, a), k=2)
    raw_resid = residual_in_span(theta, ps, hstep)
    Eps = [norm(E_map(p, x, a, hstep), hstep) for p in ps]
    E_resid = residual_in_span(theta, Eps, hstep)
    better = "E" if E_resid < raw_resid else "raw"
    print(f"{a:5.1f} {eps:9.3f} {raw_resid:11.4f} {E_resid:12.4f} {better:>8}")
