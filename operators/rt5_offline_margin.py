#!/usr/bin/env python3
"""§386: my-instrument analog of the §393 live frontier (NOT the Olympus mu^T Sig^-1 mu form —
an independent cross-check on l2_edge_sweep). Whiten the Weil form into the canonical orthonormal
basis (overlap-Gram whitening = the §391 fix that killed the frame/localization artifact), then
across sigma(=w) x density(=s) measure:
  - on-line whitened margin  = lambda_min(G^{-1/2} B_analytic G^{-1/2})   (B = pole+arch-prime, PSD=RH)
  - off-line whitened margin = lambda_min(G^{-1/2} M_beta G^{-1/2}), M_beta = the GENUINE
    functional-equation mirror-pair off-line kernel summed over band zeros:
      M_beta[j,k] = sum_{gamma in band} [ F_j(g+ib)F_k(g+ib) + F_j(g-ib)F_k(g-ib) - 2 F_j(g)F_k(g) ]
    (the exact explicit-formula off-line contribution, real part; b=beta displacement).
  - rcond(G) = conditioning of the whitening.
FORK (§393): a well-conditioned sigma where on-line >= +floor AND off-line < 0 cleanly => CHARGEABLE
spectral window (escalate). OR off-line only < 0 where rcond tiny / margin->0 as w sharpens =>
resolution/conditioning tradeoff is STRUCTURAL (spectral difficulty conservation). No-free-lunch binds.
"""
import numpy as np
import l2_edge_sweep as L

def whiten_lmin(M, G, rc=1e-12):
    # canonical orthonormal basis: G = V diag(s) V^T; keep s > rc*max; W = D^{-1/2}V^T M V D^{-1/2}
    s, V = np.linalg.eigh(G)
    keep = s > rc * s[-1]
    Vs = V[:, keep]; ss = s[keep]
    Wt = (Vs / np.sqrt(ss)).T @ M @ (Vs / np.sqrt(ss))
    return float(np.linalg.eigvalsh(Wt)[0]), int(np.sum(keep)), float(ss[0] / ss[-1])

def Fc(cs, w, r):  # F_j evaluated at (possibly complex) r
    return np.exp(-((r - cs) ** 2) / (2 * w * w)) + np.exp(-((r + cs) ** 2) / (2 * w * w))

def Mbeta(cs, w, beta):
    n = len(cs); M = np.zeros((n, n))
    zin = L.ZEROS[(L.ZEROS >= cs.min() - 2 * w) & (L.ZEROS <= cs.max() + 2 * w)]
    for g in zin:
        Fp = Fc(cs, w, g + 1j * beta); Fm = Fc(cs, w, g - 1j * beta); F0 = Fc(cs, w, g)
        M += np.real(np.outer(Fp, np.conj(Fp)) + np.outer(Fm, np.conj(Fm))) - 2.0 * np.real(np.outer(F0, F0))
    return M

BAND = (14.0, 44.0)
beta = 0.10
print(f"§386 whitened on-line vs off-line(beta={beta}) margin surface, band [14,44]. FORK per §393.")
print(f"{'w':>4} {'s':>5} {'n':>3} {'rcond':>9} {'kept':>4} {'onl_white':>11} {'offl_white':>11} {'true_lmin':>10}")
for w in (4.0, 3.0, 2.0, 1.5):
    for s in np.arange(6.0, 2.4, -0.5):
        cs = L.probes(BAND, s)
        if len(cs) < 4: continue
        f = L.build_forms(cs, w)
        G = f['overlap']; B = f['full']
        Mb = Mbeta(cs, w, beta)
        onl, kept, cond = whiten_lmin(B, G)
        offl, _, _ = whiten_lmin(Mb, G)
        tl = float(np.linalg.eigvalsh(B)[0])
        print(f"{w:4.1f} {s:5.1f} {len(cs):3d} {1/cond:9.1e} {kept:4d} {onl:+11.3e} {offl:+11.3e} {tl:+10.2e}",
              flush=True)
print("""
reading (FORK):
  CHARGEABLE window => a (w,s) with rcond NOT tiny where onl_white >= +floor AND offl_white < 0 clean.
  STRUCTURAL tradeoff => offl_white only < 0 where rcond->0 (ill-conditioned) or onl_white->0 as w sharpens.
onl_white = whitened Weil margin (on-line, PSD=RH); offl_white = whitened off-line mirror-pair kernel margin
(negative = an off-line zero drives the form negative in that basis). No-free-lunch: well-conditioned != a-priori.""")
