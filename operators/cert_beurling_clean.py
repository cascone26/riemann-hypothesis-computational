#!/usr/bin/env python3
"""§829: the CLEAN capacity-matched Beurling probe — does the Weil margin track t*=2πe^{4W}
(sampling/Beurling), or is it just Gram conditioning (the §826/§828-probe confound)?

Disambiguation by design — every row reports G-conditioning (λ_min(G), cond(G)) so the two
mechanisms can't be conflated:
  Panel 1  BEURLING SWEEP: FIXED small span (K below frame capacity), sweep 4W so t* crosses
           the window [14,44]. Prediction if Beurling: margin falls monotonically as t* rises
           across the window, WHILE G stays well-conditioned.
  Panel 2  N-CONTROL at fixed 4W (within capacity): add modes; margin should stay stable and G
           healthy → the Panel-1 collapse is from t*, not mode count.
  Panel 3  OVER-CAPACITY at the same fixed 4W: push N past frame capacity; margin collapses AND
           λ_min(G)→0 / cond→∞ → THIS collapse is the conditioning artifact (§826), visibly
           distinct from Panel 1.

Faithful float replication of operators/l2_fejer.build() with W as a parameter (self-contained to
avoid the module-state landmine that bit §806). Q = pole + arch − prime; G = ∫F_jF_k dr; margin =
min generalized eigenvalue of (Q, G). Float = REPORT-tier (measured); arb certification of the two
Panel-1 anchor points is the Handle upgrade (cert_subtstar route), run separately if the picture holds.
"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh

BAND = (14.0, 44.0)
ZEROS = np.array([
    14.134725142, 21.022039639, 25.010857580, 30.424876126, 32.935061588,
    37.586178159, 40.918719012, 43.327073281, 48.005150881, 49.773832478,
    52.970321478, 56.446247697, 59.347044003, 60.831778525, 65.112544048])

# von Mangoldt weights
def vonmangoldt(N):
    lam = np.zeros(N + 1); isc = np.zeros(N + 1, dtype=bool)
    for p in range(2, N + 1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam
NMAX = 200
LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX + 1)
LOGN = np.log(NS); LAMW = LAM[2:] / np.sqrt(NS)

RGRID = np.arange(-200.0, 200.0, 0.0025); dr = 0.0025
OMEGA = np.real(digamma(0.25 + 0.5j * RGRID)) - np.log(np.pi)

def fej(x, W):
    y = W * x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
    out[nz] = (np.sin(y[nz]) / y[nz]) ** 2
    return out
def Frow(c, r, W):
    return fej(r - c, W) + fej(r + c, W)

def build(cs, W):
    """Faithful float build: returns Q=pole+arch-prime, G=overlap. cs = centers, W = half-bandwidth."""
    n = len(cs); fourW = 4 * W
    P = np.array([Frow(c, RGRID, W) for c in cs])         # n x len(r)
    G = (P @ P.T) * dr                                     # ∫ F_jF_k dr
    arch = (P * OMEGA[None, :]) @ P.T * dr / (2 * np.pi)
    # prime side: cutoff log n < 4W (exact, finite, no tail)
    mask = LOGN < fourW
    ln = LOGN[mask]; lw = LAMW[mask]
    prime = np.zeros((n, n))
    if ln.size:
        cosu = np.cos(np.outer(ln, RGRID))                 # nprime x len(r)
        for j in range(n):
            FjP = P[j][None, :] * P                          # n x len(r)
            gp = (FjP @ cosu.T) * dr / (2 * np.pi)          # n x nprime
            prime[j, :] = 2.0 * (gp * lw[None, :]).sum(axis=1)
        prime = 0.5 * (prime + prime.T)
    # pole: F_c(i/2) = fej(i/2-c)+fej(i/2+c), fej(x)=(sin(Wx)/(Wx))^2 complex
    def fc(c):
        y = W * (0.5j - c); a = (np.sin(y) / y) ** 2
        y2 = W * (0.5j + c); b = (np.sin(y2) / y2) ** 2
        return a + b
    fp = np.array([fc(c) for c in cs])
    pole = 2.0 * np.real(np.outer(fp, fp))
    Q = pole + arch - prime
    return Q, G

def margin(cs, W):
    Q, G = build(cs, W)
    gw = np.linalg.eigvalsh(G)
    lmin_G = float(gw[0]); cond_G = float(gw[-1] / max(gw[0], 1e-300))
    try:
        m = float(eigh(Q, G, eigvals_only=True)[0])
    except np.linalg.LinAlgError:
        m = float('nan')
    return m, lmin_G, cond_G

def tstar(fourW): return 2 * np.pi * np.exp(fourW)
def zeros_below(ts): return int(((ZEROS >= 14) & (ZEROS <= 44) & (ZEROS < ts)).sum())
def zeros_in_window(): return int(((ZEROS >= 14) & (ZEROS <= 44)).sum())
def win_frac_below(ts): return float(min(max((ts - 14) / 30.0, 0.0), 1.0))
def nprime(fourW): return int((LOGN < fourW).sum())

print("=" * 100)
print("§829 CLEAN BEURLING PROBE — does the Weil margin track t*, with G-conditioning controlled?")
print("=" * 100)
Zwin = zeros_in_window()
print(f"window {BAND}, {Zwin} zeros in window; t* = 2π·e^(4W); margin = min gen-eig of (Q,G)")
print(f"CONDITIONING GUARD: λ_min(G) healthy (not →0) and cond(G) bounded ⟹ collapse is NOT the artifact\n")

# ---------- PANEL 1: Beurling sweep, FIXED K=5 span below capacity, 4W sweeps t* across window ----------
print("-" * 100)
print("PANEL 1 — BEURLING SWEEP: fixed K=5 centers, sweep 4W so t* crosses [14,44]")
print("-" * 100)
K = 5; cs1 = list(np.linspace(14, 44, K))
print(f"centers (fixed): {[round(c,2) for c in cs1]}\n")
hdr = f"{'4W':>6} {'t*':>9} {'z<t*':>6}/{Zwin} {'win%<t*':>8} {'n_prime':>8} {'margin':>13} {'λmin(G)':>11} {'cond(G)':>10}"
print(hdr)
for fourW in [1.30, 1.45, 1.60, 1.75, 1.90, 2.05]:
    W = fourW / 4; ts = tstar(fourW)
    m, lg, cg = margin(cs1, W)
    print(f"{fourW:>6.2f} {ts:>9.2f} {zeros_below(ts):>6}/{Zwin} {win_frac_below(ts):>8.2f} {nprime(fourW):>8} {m:>+13.5e} {lg:>11.3e} {cg:>10.2e}")

# ---------- PANEL 2: N-control at FIXED 4W=1.60 (t* mid-window), within capacity ----------
print("\n" + "-" * 100)
print("PANEL 2 — N-CONTROL at fixed 4W=1.60 (t*=%.1f, mid-window): add modes WITHIN capacity" % tstar(1.60))
print("prediction if Beurling (not conditioning): margin ~stable, G stays healthy" )
print("-" * 100)
W = 1.60 / 4
print(f"{'N':>4} {'centers':>28} {'margin':>13} {'λmin(G)':>11} {'cond(G)':>10}")
for N in [3, 4, 5, 6]:
    cs = list(np.linspace(14, 44, N)); m, lg, cg = margin(cs, W)
    print(f"{N:>4} {str([round(c,1) for c in cs]):>28} {m:>+13.5e} {lg:>11.3e} {cg:>10.2e}")

# ---------- PANEL 3: OVER-CAPACITY at same fixed 4W=1.60 — expose the conditioning artifact ----------
print("\n" + "-" * 100)
print("PANEL 3 — OVER-CAPACITY at fixed 4W=1.60: push N past frame capacity")
print("prediction: margin collapses AND λmin(G)→0 / cond(G)→∞ ⟹ THIS collapse = conditioning artifact (§826)")
print("-" * 100)
print(f"{'N':>4} {'margin':>13} {'λmin(G)':>11} {'cond(G)':>10}  {'verdict'}")
for N in [7, 9, 12, 16]:
    cs = list(np.linspace(14, 44, N)); m, lg, cg = margin(cs, W)
    verdict = "ARTIFACT (G ill-conditioned)" if (lg < 1e-6 or cg > 1e8) else "healthy"
    print(f"{N:>4} {m:>+13.5e} {lg:>11.3e} {cg:>10.2e}  {verdict}")

print("\n" + "=" * 100)
print("READ (main-session interpretation added after run):")
print(" - If Panel 1 margin falls as t* rises WHILE λmin(G) stays healthy → Beurling mechanism SHOWN (real).")
print(" - If Panel 2 margin is stable under added modes → the collapse is NOT mode-count/conditioning.")
print(" - If Panel 3 margin collapse coincides with λmin(G)→0 → that (and only that) is the §826 artifact.")
print("=" * 100)
