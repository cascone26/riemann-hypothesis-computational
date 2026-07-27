#!/usr/bin/env python3
"""§365: THE EXCLUSION CERTIFICATE, float proof-of-concept (arb hardening after lit-gate).

Two-sided test at the Nyquist/rank row of a window:
  UPPER (zero-free, primes+arch+pole only): u = analytic-side near-null eigenvector of B_full;
    U := <u, B_full u> computed from pole+arch-prime (NO zeros). Small U => the analytic form is
    near-singular in direction u.
  LOWER (the fake worlds): for candidate off-line displacement delta, L(delta) := min over unit
    vectors of <v, B_zero^{fake}(delta) v> at the SAME probes = smallest eigenvalue the zero side
    is FORCED to have if a zero sits at 1/2+delta+i*gamma in-window. (Here measured directly.)
  EXCLUSION LOGIC: B_full == B_zero (explicit formula, residual ~1e-10). So if the TRUE world had
    an off-line zero at displacement >= delta in-window, B_full would inherit L(delta) > 0 in every
    direction -> in particular <u,B_full u> >= L(delta). Thus  U < L(delta)  is INCONSISTENT with
    such a zero: measuring U (primes only) < c(delta) EXCLUDES off-line zeros of displacement >=
    delta effectively in-window. This script measures U and L(delta) and prints the excluded delta.
Caveat tracked: localization — B_zero sums ALL zeros; the fake world moves ONE in-window zero and
keeps the rest on-line. Out-of-window leakage = the honest gap (flagged, quantified by the residual
of the on-line-only reconstruction vs full).
"""
import numpy as np
import l2_edge_sweep as L

def zside_complex(cs, w, zlist):
    n = len(cs); B = np.zeros((n, n))
    for z in zlist:
        Fv = np.exp(-((z - cs) ** 2) / (2 * w * w)) + np.exp(-((z + cs) ** 2) / (2 * w * w))
        B += 2.0 * np.real(np.outer(Fv, np.conj(Fv)))
    return B

print("EXCLUSION CERTIFICATE (float PoC): U = <u,B_full u> (primes only) vs L(delta) (fake worlds)")
print(f"{'w':>4} {'s':>5} {'U (primes)':>11} {'true lmin':>10} | " +
      " | ".join(f"L({d})" for d in (0.05, 0.1, 0.2)) + "  excl_delta")
for w in (3.0, 2.0, 1.5):
    # rank row: first s with true zero-side lmin < 1e-12
    srow = None
    for s in np.arange(6.0, 1.19, -0.1):
        cs = L.probes(L.BAND, s)
        if len(cs) < 4: continue
        f = L.build_forms(cs, w)
        if float(np.linalg.eigvalsh(f['zero'])[0]) < 1e-12:
            srow = s; break
    cs = L.probes(L.BAND, srow)
    f = L.build_forms(cs, w)
    B = f['full']                        # = pole+arch-prime, zero-free
    lam, U_ = np.linalg.eigh(B)
    u = U_[:, 0]                         # analytic near-null direction
    Uval = float(u @ B @ u)              # primes-only scalar (== lam[0])
    true_lmin = float(np.linalg.eigvalsh(f['zero'])[0])
    Ls = []
    for d in (0.05, 0.1, 0.2):
        # fake world: move EACH in-window on-line zero to 1/2+d (test the strongest: the resonant one)
        # here: displace the single zero whose F-projection on u is largest (the detected one)
        gproj = [(abs(float((np.exp(-((g - cs)**2)/(2*w*w)) + np.exp(-((g + cs)**2)/(2*w*w))) @ u)), g)
                 for g in L.ZEROS]
        gsel = max(gproj)[1]
        fake = [complex(g, 0.0) for g in L.ZEROS if g != gsel] + [complex(gsel, d)]
        Lf = float(np.linalg.eigvalsh(zside_complex(cs, w, fake))[0])
        Ls.append(Lf)
    # excluded delta: smallest d with L(d) > U (U below the forced floor => that zero can't exist)
    excl = None
    for d, Lf in zip((0.05, 0.1, 0.2), Ls):
        if Lf > max(Uval, 0) + 1e-15:
            excl = d; break
    print(f"{w:4.1f} {srow:5.1f} {Uval:+11.3e} {true_lmin:+10.2e} | " +
          " | ".join(f"{Lf:8.2e}" for Lf in Ls) + f"   >= {excl}", flush=True)
print("""
Reading: U is computed from PRIMES ONLY (analytic side). Where U < L(delta), a true off-line zero at
displacement delta in-window is EXCLUDED (it would force lmin >= L(delta) > U, contradicting the
measured primes-only U). This is the zero-free half; the rigorous version needs (a) arb interval U
(upper), (b) an analytic lower bound on L(delta) [the c(delta)=kappa*delta^2 law, transversality of
the Im-mode], (c) out-of-window leakage charge. Lit-gate pending before any novelty claim.""")
