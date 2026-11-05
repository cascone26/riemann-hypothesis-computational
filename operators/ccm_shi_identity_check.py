#!/usr/bin/env python3
"""HP RH center, cycle 1 (2026-11-04): is CCM's Weil matrix Q (arXiv:2511.22755, Lemma 5.1 form,
as built in operators/ccm_exact_reproduce.py) EXACTLY Shi's S_{N,L} (arXiv:2609.04908 Eqs 3.35-3.39)?
Independent oracle: the zero side.  Shi defines S_mn = Z(G_mn) = sum over nontrivial zeros of the
Fourier transform of the real-polarized autocorrelation of f_n = L^{-1/2} e^{i nu_n x} 1_[0,L],
nu_n = 2 pi n / L.  Under RH-consistent data (all 100k zeros on the line) that sum is
   Z_mn = sum_k Re[ fh_m(g) conj fh_n(g) + fh_{-m}(g) conj fh_{-n}(g) ] at g = gamma_k,
   fh_n(t) = L^{-1/2} (e^{i t L} - 1) / (i (nu_n + t)).
Tail beyond gamma_max corrected with the Riemann-von Mangoldt density (1/2pi) log(t/2pi).
"""
import sys, numpy as np, mpmath as mp
sys.path.insert(0, 'operators')
from ccm_exact_reproduce import build

lam2 = float(sys.argv[1]) if len(sys.argv) > 1 else 14.0
N = int(sys.argv[2]) if len(sys.argv) > 2 else 4
L = np.log(lam2)
idx = np.arange(-N, N + 1)
nu = 2 * np.pi * idx / L

g = np.loadtxt('data/zeros_100k.txt', usecols=-1) if True else None
g = np.asarray(g, float)
def fh(t):  # rows: t, cols: n
    t = t[:, None]
    return (np.exp(1j * t * L) - 1) / (1j * (nu[None, :] + t)) / np.sqrt(L)
Z = np.zeros((len(idx), len(idx)))
for chunk in np.array_split(g, 50):
    F = fh(chunk)
    Fm = F[:, ::-1]  # fh_{-n}
    Z += np.real(F.T @ F.conj() + Fm.T @ Fm.conj())
# tail: average of Re[...] over the phase e^{itL} -> |e^{itL}-1|^2 averages to 2; entries ~ 2*2/(L t^2)
T = g[-1]
tail = mp.quad(lambda t: mp.log(t / (2 * mp.pi)) / (2 * mp.pi) * 4 / (L * t ** 2), [T, mp.inf])
Zt = Z + float(tail)  # leading-order tail (same for all entries to O(1/T^2))

Q, Lm, _ = build(mp.sqrt(lam2), N)
Qn = np.array([[float(Q[i, j]) for j in range(len(idx))] for i in range(len(idx))])

np.set_printoptions(precision=6, suppress=True, linewidth=160)
print(f"lambda^2={lam2} L={L:.6f} N={N}  zeros used={len(g)} gamma_max={T:.1f} tail={float(tail):.3e}")
print("zero-side Z (corner 4x4):\n", Zt[:4, :4])
print("CCM Q (corner 4x4):\n", Qn[:4, :4])
R = Qn / Zt
print("ratio Q/Z  min/max/mean:", R.min(), R.max(), R.mean())
c = (Qn * Zt).sum() / (Zt * Zt).sum()
print("best scalar c (Q ~ c Z):", c, "  max|Q - cZ| =", np.abs(Qn - c * Zt).max(), "  max|Q|=", np.abs(Qn).max())
print("max|Q - Z| =", np.abs(Qn - Zt).max())
np.save(f'results/hp_center/ccm_shi_Q_N{N}.npy', Qn); np.save(f'results/hp_center/ccm_shi_Z_N{N}.npy', Zt)
