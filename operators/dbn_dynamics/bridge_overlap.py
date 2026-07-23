#!/usr/bin/env python3
"""§178: k_λ BRIDGE increment 1 — overlap(k_λ, ψ_CvS) at CvS/mpmath precision.
The §176 cascade lives on the CvS numerical near-null eigenvector ψ. Connes' analytical k_λ=E(h_λ) (prolate combo) is the
object with a rigorous handle (Fuchs). §177 feasibility: in the OLD float64 Suzuki operator overlap(v,k_λ)=1.0000 at λ=1.8
but that operator saturates ~1e-7 so it's blind past λ≈1.8. Here: compute overlap in the CvS mpmath operator (unsaturated),
extending the confirmation to larger λ (c=13 ⇒ λ=3.6). If overlap→1, k_λ IS the CvS near-null vector ⇒ the cascade transfers.

ψ_CvS even-sector coeff vector c_k ⇒ function f_ψ(t)=c_0 + Σ_{k≥1}√2 c_k cos(πk t/L), t∈[−L,L], L=log c.
k_λ(s): analytical prolate vector on s∈[−log λ,log λ]=[−L/2,L/2] (λ=√c). Convention map tested: s=t/2 (CvS t=2 log u, c=λ²).
Anchor: at small c overlap must reproduce stage3's ~1.0000 or the convention is wrong. k_λ from float64 prolates (~1e-15,
plenty for an overlap-→1 test; still 8 orders past the float64 OPERATOR floor 1e-7)."""
import sys, time
import numpy as np
import mpmath as mp
sys.path.insert(0, ".venv_cvs/lib/python3.14/site-packages")
import connes_cvs.operator as O
from numpy.polynomial.legendre import legval

# ---------------- analytical k_λ (numpy float64, from stage3_theta_vs_k.py) ----------------
def _y2(Nleg):
    Y = np.zeros((Nleg + 2, Nleg + 2))
    for n in range(Nleg + 1):
        if n + 1 <= Nleg + 1: Y[n + 1, n] = (n + 1) / (2 * n + 1) * np.sqrt((2 * n + 1) / (2 * n + 3))
        if n - 1 >= 0:        Y[n - 1, n] = n / (2 * n + 1) * np.sqrt((2 * n + 1) / (2 * n - 1))
    return (Y[:Nleg, :Nleg]) @ (Y[:Nleg, :Nleg])

def k_lambda_on(s_grid, lmbda):
    Nleg = max(200, int(6 * 2 * np.pi * lmbda ** 2)); cband = 2 * np.pi * lmbda ** 2; n = np.arange(Nleg)
    A = np.diag((n * (n + 1)).astype(float)) + cband ** 2 * _y2(Nleg); A = (A + A.T) / 2
    w, V = np.linalg.eigh(A)
    par = np.array([(-1) ** int(np.argmax(np.abs(V[:, k]))) for k in range(V.shape[1])])
    ei = np.where(par > 0)[0]; i0, i4 = ei[0], ei[2]
    a0 = 1.0; b0 = (-V[0, i0] / V[0, i4] if abs(V[0, i4]) > 1e-30 else 0.0)
    def coeff(col): return np.array([col[m] * np.sqrt((2 * m + 1) / 2) for m in range(Nleg)])
    c0, c4 = coeff(V[:, i0]), coeff(V[:, i4])
    def h_lam(xv):
        y = xv / lmbda; out = np.zeros_like(xv, float); m = np.abs(y) <= 1
        out[m] = a0 * legval(y[m], c0) + b0 * legval(y[m], c4); return out
    u = np.exp(s_grid); nmax = int(np.floor(lmbda ** 2)) + 1
    ssum = np.zeros_like(u)
    for nn in range(1, nmax): ssum += h_lam(nn * u)
    return np.sqrt(u) * ssum

# ---------------- CvS near-null eigenvector ψ (even-sector coeffs) ----------------
def cvs_ground(c, N, T, dps):
    mp.mp.dps = dps
    if O.HAS_FLINT:
        from flint import ctx as _fctx; _fctx.prec = int(dps * 3.5)
    L = mp.log(c); pdata, _ = O.prime_powers_up_to(int(c))
    idx = list(range(-N, N + 1)); psiv = {}; psid = {}
    for m in idx:
        x = mp.mpf(m); O._hplus_cache_clear()
        psiv[m] = O.psi_prime(x, L, pdata) + O.psi_pole(x, L) + O.psi_arch(x, L, T, dps)
        psid[m] = O.psi_prime_deriv(x, L, pdata) + O.psi_pole_deriv(x, L) + O.psi_arch_deriv(x, L, T, dps)
        O._hplus_cache_clear()
    DIM = 2 * N + 1; Q = mp.matrix(DIM, DIM)
    for i in range(DIM):
        m = i - N
        for j in range(DIM):
            nn = j - N
            Q[i, j] = psid[nn] if m == nn else (psiv[m] - psiv[nn]) / (m - nn)
    for i in range(DIM):
        for j in range(i + 1, DIM):
            a = (Q[i, j] + Q[j, i]) / 2; Q[i, j] = a; Q[j, i] = a
    Vv = mp.matrix(DIM, N + 1); Vv[N, 0] = mp.mpf(1); r2 = 1 / mp.sqrt(2)
    for k in range(1, N + 1): Vv[N + k, k] = r2; Vv[N - k, k] = r2
    Qe = Vv.T * Q * Vv
    eigs, vecs = mp.eigsy(Qe)
    k0 = min(range(N + 1), key=lambda i: eigs[i])
    ev = [float(vecs[r, k0]) for r in range(N + 1)]       # even-sector coeffs c_0..c_N (float is plenty for overlap)
    return float(mp.log10(abs(eigs[k0]))), ev, float(L)

def f_psi(x, ev, L):
    """CvS eigenFUNCTION on x∈[0,L] (paper Eq5: basis U_n=L^{-1/2}exp(2πinx/L), period L). Even sector (η_{-k}=η_k):
    η(x) = c_0 + Σ_{k≥1} √2 c_k cos(2π k x / L)."""
    out = np.full_like(x, ev[0], float)
    for k in range(1, len(ev)):
        out += np.sqrt(2) * ev[k] * np.cos(2 * np.pi * k * x / L)
    return out

def overlap(A, B, t):
    A = A / np.sqrt(np.trapezoid(A * A, t)); B = B / np.sqrt(np.trapezoid(B * B, t))
    ov = np.trapezoid(A * B, t)
    if ov < 0: B = -B; ov = -ov
    return ov, np.sqrt(np.trapezoid((A - B) ** 2, t))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    dps = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    T = int(sys.argv[3]) if len(sys.argv) > 3 else 400
    clist = [int(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4 else [3, 5, 8, 13]
    print(f"§178 k_λ↔ψ_CvS overlap  N={N} dps={dps} T={T}  (paper convention: x∈[0,L] period L, k_λ shifted s=x−L/2)\n")
    print(f"{'c':>4} {'λ':>7} {'log10λmin':>10} {'overlap':>10} {'dist':>9}   {'best-α scan':>18}")
    for c in clist:
        t0 = time.time()
        lmin_log, ev, L = cvs_ground(c, N, T, dps)
        lam = np.sqrt(c)
        M = 6000
        x = np.linspace(0.0, L, M + 1)          # CvS eigenfunction lives on [0,L]
        fp = f_psi(x, ev, L)
        s = x - L / 2.0                          # Connes variable: s=log u ∈ [−L/2,L/2] = [−log λ, log λ]
        kA = k_lambda_on(s, lam)
        ov, d = overlap(fp, kA, x)
        # convention robustness: scan a residual affine shift/scale s = α(x−L/2), find best α (should peak ≈ 1.0)
        best = (0.0, 0.0)
        for alpha in np.linspace(0.60, 1.60, 51):
            kk = k_lambda_on(np.clip(alpha * s, -np.log(lam), np.log(lam)), lam)
            ova, _ = overlap(fp, kk, x)
            if ova > best[1]: best = (alpha, ova)
        print(f"{c:>4} {lam:>7.3f} {lmin_log:>10.2f} {ov:>10.6f} {d:>9.4f}   α*={best[0]:.3f} ov*={best[1]:.6f}  ({time.time()-t0:.0f}s)", flush=True)
    print("\nREAD: overlap→1 (α*≈1.0) confirms k_λ = the CvS near-null eigenvector at mpmath precision, EXTENDING §177's")
    print("float64 overlap=1.0000 (only reachable λ≤1.8) up to λ=√13≈3.6 where the old operator was saturated/blind.")
    print("If so, the §176 cascade lives on the ANALYTICAL k_λ ⇒ description → proof-ingredient (ties to rigorous Fuchs).")
