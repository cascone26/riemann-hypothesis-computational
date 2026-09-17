"""
DECISIVE step-(b) test: does the finite Weil minimizer θ_x (lowest eigenvector of the screw-g form on
[-a,a], additive log-variable x) match Connes' EXACT limiting object k(u)=E(h)(u) under u=e^x, and does
the overlap GROW toward 1 as the window a grows?  Baseline: raw ground prolate ps_0 (known to NOT capture
θ_x per step_b_prolate_span).  If |<θ, k(e^x)>| → 1 while |<θ, ps_0>| does not, Connes' k_λ mechanism is
confirmed on our own numerics: θ_x → k_λ, the object whose Mellin transform is Ξ.
"""
import numpy as np
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def h(u):
    return (np.pi / 2) * u**2 * (2 * np.pi * u**2 - 3) * np.exp(-np.pi * u**2)

def k_connes(u, Nmax=4000):
    u = np.atleast_1d(u).astype(float)
    s = np.zeros_like(u)
    for n in range(1, Nmax + 1):
        term = h(n * u)
        s += term
        if np.max(np.abs(term)) < 1e-18:
            break
    return np.sqrt(u) * s

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

def ground_prolate(a, x, hstep, c_bw):
    D = x[:, None] - x[None, :]
    K = np.where(np.abs(D) < 1e-12, c_bw / np.pi, np.sin(c_bw * D) / (np.pi * D)) * hstep
    ev, V = np.linalg.eigh(K)
    g = V[:, -1]
    if g[np.argmax(np.abs(g))] < 0: g = -g
    return g / np.sqrt(np.sum(g * g) * hstep)

def norm(f, hstep):
    return f / np.sqrt(np.sum(f * f) * hstep)

def overlap(f, g, hstep):
    return abs(np.sum(f * g) * hstep)

print("θ_x (Weil minimizer) vs Connes k(e^x)  and vs raw ground prolate ps_0, overlap vs window a")
print(f"{'a':>5} {'ε(min-eig)':>11} {'|<θ,k(e^x)>|':>13} {'|<θ,ps0>|':>11} {'winner':>8}")
for a in [1.0, 1.5, 2.0, 2.5, 3.0]:
    M = min(1400, max(600, int(700 / a)))
    x, hstep, theta, eps = weil_min(a, M)
    # Connes k in the log variable u = e^x ; k as an L^2(dx) function
    kx = norm(k_connes(np.exp(x)), hstep)
    # even part of theta (k(e^x) is not symmetric in x; compare even components on equal footing)
    ps0 = ground_prolate(a, x, hstep, c_bw=max(1.0, a))
    ov_k = overlap(theta, kx, hstep)
    ov_p = overlap(theta, ps0, hstep)
    win = "k" if ov_k > ov_p else "ps0"
    print(f"{a:5.1f} {eps:11.3e} {ov_k:13.4f} {ov_p:11.4f} {win:>8}")
