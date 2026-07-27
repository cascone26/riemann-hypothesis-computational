#!/usr/bin/env python3
"""§364: THE OFF-LINE-ZERO SIGNATURE, MEASURED (the §363 gate clause with numbers).

Move the first zero rho_1 = 1/2 + i*14.1347 off the line to beta = 1/2 + delta (quadruple
beta+i*gamma, 1-beta+i*gamma and conjugates; for our even real probes the zero-side sum
becomes 2*Re[F_j(gamma+i*delta)F_k(gamma-i*delta)] + mirror). Measure:
  D = B_zero_fake - B_zero_true  (the non-L2 mode's imprint on the form)
vs the true margin, per (w, s): the DETECTION FRONTIER = deepest s at which ||D|| > margin
(where a certificate of the true form would CONTRADICT the fake world — i.e., where our
existing theorems already exclude such a zero's local effect... stated carefully: where the
fake and true worlds differ by more than the certified margin).
Also: D's minimizer projection vs the H-L variance band (the S-rigidity noise floor §352C):
the scale where the violation exceeds arithmetic noise = where the corner inequality would bite.
Exploratory/quantification — no pre-registered numeric (the OBJECT is the deliverable).
"""
import numpy as np
import l2_edge_sweep as L

def zero_side(cs, w, zeros_complex):
    """B[j,k] = 2 Re sum_z F_j(z) conj(F_k(z)) for complex zero ordinates z = gamma + i*delta
    (with the beta and 1-beta pair both included via the +/- delta list)."""
    n = len(cs)
    B = np.zeros((n, n))
    for z in zeros_complex:
        Fv = np.exp(-((z - cs) ** 2) / (2 * w * w)) + np.exp(-((z + cs) ** 2) / (2 * w * w))
        B += 2.0 * np.real(np.outer(Fv, np.conj(Fv)))
    return B

print("OFF-LINE-ZERO SIGNATURE: first zero moved to beta = 1/2 + delta")
print(f"{'delta':>6} {'w':>4} {'s':>5} {'true margin':>11} {'||D||':>10} {'D/margin':>9} {'detectable':>10}")
for delta in (0.05, 0.1, 0.2):
    for w in (2.0, 1.5):
        frontier = None
        for s in np.arange(6.0, 1.19, -0.1):
            cs = L.probes(L.BAND, s)
            if len(cs) < 4: continue
            f = L.build_forms(cs, w)
            tm = float(np.linalg.eigvalsh(f['full'])[0])
            if tm < 1e-13: break
            true_z = [complex(g, 0.0) for g in L.ZEROS]
            fake_z = ([complex(L.ZEROS[0], delta), complex(L.ZEROS[0], -delta)]
                      + [complex(g, 0.0) for g in L.ZEROS[1:]])
            # true first zero contributes once; fake contributes the beta/(1-beta) pair (each once)
            Bt = zero_side(cs, w, true_z)
            Bf = zero_side(cs, w, fake_z) - zero_side(cs, w, [complex(L.ZEROS[0], 0.0)])
            # careful bookkeeping: fake replaces ONE on-line zero by the off-line PAIR minus... the
            # explicit formula counts each zero once: fake world = pair {beta, 1-beta} at gamma_1
            # = the two complex ordinates above; true world = one ordinate at gamma_1 (mult 1)...
            # but a beta/(1-beta) pair has TWO zeros where the true world has TWO on-line zeros only
            # if gamma_1 had multiplicity 2. Simplest honest comparison: ONE simple zero moved to
            # beta (delta-ordinate) vs staying at 1/2: D below uses the single-ordinate difference.
            Fv_t = np.exp(-((L.ZEROS[0] - cs) ** 2) / (2 * w * w)) + np.exp(-((L.ZEROS[0] + cs) ** 2) / (2 * w * w))
            zc = complex(L.ZEROS[0], delta)
            Fv_f = np.exp(-((zc - cs) ** 2) / (2 * w * w)) + np.exp(-((zc + cs) ** 2) / (2 * w * w))
            D = 2.0 * (np.real(np.outer(Fv_f, np.conj(Fv_f))) - np.real(np.outer(Fv_t, Fv_t)))
            nD = float(np.linalg.norm(D, 2))
            det = nD > tm
            if det and frontier is None:
                frontier = s
            if s in (5.0, 4.5, 4.0) or (det and frontier == s):
                print(f"{delta:6.2f} {w:4.1f} {s:5.1f} {tm:11.3e} {nD:10.3e} {nD/tm:9.2e} "
                      f"{'YES' if det else 'no':>10}", flush=True)
        print(f"   -> delta={delta} w={w}: detection frontier s = {frontier} "
              f"(rank-limit region begins ~{s:.1f})", flush=True)
