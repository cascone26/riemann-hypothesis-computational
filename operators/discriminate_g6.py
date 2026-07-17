#!/usr/bin/env python3
"""
discriminate_g6.py — N=320,448 test to distinguish γ₆ j=29(127) vs j=30(151).

Both are stable through N=256:
  j=29 (127, prime): N64=-11.16, N128=-11.00, N192=-10.85, N256=-10.54  Δ=+0.62 (slow rise)
  j=30 (151, prime): N64=-11.10, N128=-11.32, N192=-10.69, N256=-11.54  Δ=-0.44 (oscillating)

Prediction:
  j=30 (151) = TRUE RESONANCE → oscillates around -11 through N=320,448
  j=29 (127) = PRIME STABILITY FLOOR → slowly diverges; log approaches -9 by N=512

Verification: if j=29 log(N=320) > -9.5 OR continues monotonic rise, it's diverging.
If j=30 log(N=320) < -10, true resonance confirmed.

Also test γ₇ j=37(293) at N=320,448 for same confirmation.
Also test γ₉ j=44(317) at N=320 (was probable — verify).
Also test γ₄ j=30(487) at N=448 — does 487 diverge, leaving 491 as unique?
"""
import math, numpy as np, scipy.linalg
from scipy.special import digamma

TWO_PI = 2 * math.pi
LOGPI = math.log(math.pi)
G_PER_N = 900


def make_sieve(MAX_LSQ):
    sv = np.ones(MAX_LSQ + 1, dtype=bool)
    sv[0] = sv[1] = False
    for i in range(2, int(MAX_LSQ ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    return np.where(sv)[0]


sv5000_sieve = None
primes5000 = make_sieve(5000)


def compute_ev0(lsq, N):
    G = G_PER_N * N
    L = 2 * math.log(math.sqrt(lsq))
    ns = np.arange(-N, N + 1)
    D = len(ns)
    R = TWO_PI * N / L + 50.0
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]
    w = np.full(G, h)
    w[0] = w[-1] = 0.5 * h
    Om = np.real(digamma(0.25 + 0.5j * grid)) - LOGPI
    c = TWO_PI * ns[:, None] / L
    d = c - grid[None, :]
    sm = np.abs(d) < 1e-12
    E = np.where(sm, np.sqrt(L) + 0j,
                 (1 / np.sqrt(L)) * np.exp(0.5j * d * L) * 2 * np.sin(0.5 * d * L) / d)
    arch = (E * (w * Om) @ E.conj().T).real / TWO_PI
    sh = np.sinh(L / 4) ** 2
    ni = ns[:, None]
    nj = ns[None, :]
    B = 32 * L * sh * (L ** 2 - 16 * math.pi ** 2 * ni * nj) / (
        (L ** 2 + 16 * math.pi ** 2 * ni ** 2) * (L ** 2 + 16 * math.pi ** 2 * nj ** 2))
    p_in = primes5000[primes5000 <= lsq]
    ll = []
    ww = []
    for p in p_in:
        lp = math.log(p)
        pk = p
        while pk <= lsq:
            ll.append(math.log(pk))
            ww.append(lp / math.sqrt(pk))
            pk *= p
    if not ll:
        return 0.0
    lkpp = np.array(ll)
    wts_p = np.array(ww)
    K = len(lkpp)
    S = np.zeros(D)
    Sd = np.zeros(D)
    for k0 in range(0, K, 10000):
        k1 = min(k0 + 10000, K)
        ph = TWO_PI * np.outer(ns, lkpp[k0:k1]) / L
        S += np.sin(ph) @ wts_p[k0:k1]
        Sd += np.cos(ph) @ (2 * wts_p[k0:k1] * (1 - lkpp[k0:k1] / L))
    dn = ns[:, None] - ns[None, :]
    Si = S[:, None]
    Sj = S[None, :]
    with np.errstate(divide='ignore', invalid='ignore'):
        P = np.where(dn != 0, (Sj - Si) / (math.pi * dn), 0.0)
    np.fill_diagonal(P, Sd)
    tau = (arch + B - P)
    tau = (tau + tau.T) / 2
    return float(scipy.linalg.eigvalsh(tau, subset_by_index=[0, 0])[0])


def test_series(lsq, Nvals, label):
    print(f"\n--- {label} (lsq={lsq}) ---")
    row = {}
    for N in Nvals:
        e = compute_ev0(lsq, N)
        log_e = math.log10(abs(e)) if abs(e) > 0 else -99
        row[N] = log_e
        print(f"  N={N:>4}: {log_e:>8.3f}", flush=True)
    Nlist = sorted(Nvals)
    delta = row[Nlist[-1]] - row[Nlist[0]]
    stable = abs(delta) < 2.5 and row[Nlist[-1]] < -10
    verdict = "STABLE (TRUE RESONANCE)" if stable else f"DIVERGING (Δ={delta:+.2f})"
    print(f"  → {verdict}")
    return row, stable


# ====================================================================
# γ₆ DISCRIMINATION: j=29(127) vs j=30(151) at N=64..448
# Includes the N=64..256 we know, then adds N=320 and N=448
# ====================================================================
GAMMA6 = 37.586178158651
print("=" * 65)
print(f"γ₆ DISCRIMINATION: j=29(127) vs j=30(151)")
print(f"Known through N=256:")
print(f"  j=29(127): N64=-11.16, N128=-11.00, N192=-10.85, N256=-10.54  Δ=+0.62")
print(f"  j=30(151): N64=-11.10, N128=-11.32, N192=-10.69, N256=-11.54  Δ=-0.44")
print(f"Now extending to N=320, N=448...")
print("=" * 65)

# Test at N=320 and N=448 only (we already have N=64..256)
r127, s127 = test_series(127, [320, 448], "γ₆ j=29 (lsq=127, prime, slow-rise?)")
r151, s151 = test_series(151, [320, 448], "γ₆ j=30 (lsq=151, prime, true-resonance?)")

print("\n--- γ₆ VERDICT ---")
print(f"j=29(127): N320={r127.get(320,'?'):.3f}, N448={r127.get(448,'?'):.3f}")
print(f"j=30(151): N320={r151.get(320,'?'):.3f}, N448={r151.get(448,'?'):.3f}")
if r151.get(320, 0) < -10 and r127.get(320, 0) > -9:
    print("CLEAR DISCRIMINATION: j=30(151) is TRUE RESONANCE, j=29(127) is prime floor.")
elif r151.get(320, 0) < -10 and r127.get(320, 0) < -10:
    print("BOTH STILL STABLE at N=320 — need N=512+ or different discriminator")
else:
    print("Check trajectories manually.")

# ====================================================================
# γ₇ CONFIRMATION: j=37(293) at N=320, N=448
# Prior: N64=-11.89, N128=-11.97, N192=-10.82, N256=-11.45, Δ=-0.44
# ====================================================================
GAMMA7 = 40.918719012680879
print("\n" + "=" * 65)
print(f"γ₇ CONFIRMATION: j=37 (lsq=293)")
print(f"Prior: N64=-11.89, N128=-11.97, N192=-10.82, N256=-11.45")
r293, s293 = test_series(293, [320, 448], "γ₇ j=37 (lsq=293, prime)")
print(f"\nγ₇ j=37(293) extended: {'CONFIRMED STABLE' if s293 else 'DIVERGING'}")

# ====================================================================
# γ₉ CONFIRMATION: j=44(317) at N=320
# Prior: was "probable" — need N=320 verification
# ====================================================================
GAMMA9 = 48.005150881167
j44_g9 = round(math.exp(TWO_PI * 44 / GAMMA9))
print("\n" + "=" * 65)
print(f"γ₉ CONFIRMATION: j=44 → lsq={j44_g9}")
r9, s9 = test_series(j44_g9, [64, 128, 192, 256, 320], "γ₉ j=44")
print(f"\nγ₉ j=44(lsq={j44_g9}): {'CONFIRMED' if s9 else 'DIVERGING'}")

# ====================================================================
# γ₄ COMPETITOR: lsq=487 at N=448
# Prior: BOTH 487 and 491 stable through N=320.
# Does 487 eventually diverge?
# ====================================================================
print("\n" + "=" * 65)
print(f"γ₄ COMPETITOR lsq=487 at N=448")
print(f"Prior: both 487 and 491 stable through N=320")
r487, s487 = test_series(487, [448], "γ₄ competitor lsq=487")
r491, s491 = test_series(491, [448], "γ₄ true lsq=491")
print(f"\nγ₄ at N=448: lsq=487 → {'stable' if s487 else 'DIV'}, lsq=491 → {'stable' if s491 else 'DIV'}")

print("\n\n=== FINAL STATUS ===")
print(f"γ₄ (491): confirmed={s491}, competitor(487): {not s487}")
print(f"γ₆ (151): confirmed={s151}, floor(127): {not s127}")
print(f"γ₇ (293): confirmed={s293}")
print(f"γ₉ ({j44_g9}): confirmed={s9}")
