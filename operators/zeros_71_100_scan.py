#!/usr/bin/env python3
"""
zeros_71_100_scan.py

CvS spectral operator scan for Riemann zeros γ₇₁–γ₁₀₀.
Continuation of zeros_61_100_scan.py (which completed γ₆₁–γ₇₀).
Zero values from mpmath zetazero() at 25 decimal places.
"""
import math, numpy as np, scipy.linalg, sys
from scipy.special import digamma

TWO_PI = 2 * math.pi
LOGPI = math.log(math.pi)
G_PER_N = 900

ZEROS_71_100 = [
    (71, 184.874467848387),
    (72, 185.598783677707),
    (73, 187.228922583502),
    (74, 189.416158656017),
    (75, 192.026656360714),
    (76, 193.079726603846),
    (77, 195.265396679529),
    (78, 196.876481840958),
    (79, 198.015309676252),
    (80, 201.264751943704),
    (81, 202.493594514141),
    (82, 204.189671803105),
    (83, 205.394697202163),
    (84, 207.906258887806),
    (85, 209.576509716856),
    (86, 211.690862595365),
    (87, 213.347919359713),
    (88, 214.547044783491),
    (89, 216.169538508264),
    (90, 219.067596349021),
    (91, 220.714918839314),
    (92, 221.430705554693),
    (93, 224.007000254604),
    (94, 224.983324669582),
    (95, 227.421444279679),
    (96, 229.337413305525),
    (97, 231.250188700499),
    (98, 231.987235253180),
    (99, 233.693404178908),
    (100, 236.524229665816),
]


def make_sieve(MAX_LSQ):
    sv = np.ones(MAX_LSQ + 1, dtype=bool)
    sv[0] = sv[1] = False
    for i in range(2, int(MAX_LSQ ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    return sv, np.where(sv)[0]


sv30000, primes30000 = make_sieve(30000)


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
    p_in = primes30000[primes30000 <= lsq]
    ll = []
    ww = []
    for p in p_in:
        lp = math.log(p)
        pk = p
        while pk <= lsq:
            ll.append(math.log(pk))
            ww.append(lp / math.sqrt(pk))
            pk *= p
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
    ev0 = float(scipy.linalg.eigvalsh(tau, subset_by_index=[0, 0])[0])
    return ev0


def find_resonant_prime(gamma, p_min=50, p_max=20000, top_k=8):
    candidates = []
    for p in primes30000[(primes30000 >= p_min) & (primes30000 <= p_max)]:
        L = math.log(p)
        j = round(gamma * L / TWO_PI)
        if j < 1:
            continue
        gamma_eff = TWO_PI * j / L
        err = abs(gamma_eff - gamma)
        candidates.append((err, p, j, gamma_eff))
    candidates.sort()
    return candidates[:top_k]


def stability_test(p, N_list):
    results = []
    prev = None
    for N in N_list:
        ev0 = compute_ev0(p, N)
        log_ev0 = math.log10(abs(ev0)) if abs(ev0) > 1e-99 else -99.0
        delta = (log_ev0 - prev) if prev is not None else None
        results.append((N, log_ev0, delta))
        prev = log_ev0
    return results


def classify_stability(results, p, max_p=5000):
    last_n, last_log, last_delta = results[-1]
    if last_log >= -10:
        return "div"
    if last_delta is not None and last_delta > 2.5:
        return "div"
    if len(results) >= 2 and results[-2][1] < -10 and last_log < -10:
        return "STABLE***"
    return "STABLE(1side,Δ={:+.2f})".format(last_delta if last_delta is not None else 0)


N_SEQ_LARGE = [64, 128, 192, 256]      # p >= 5000
N_SEQ_SMALL = [64, 128, 192, 320]       # p < 5000

confirmed = []
failed = []

print(f"CvS Operator scan: γ₇₁–γ₁₀₀ (continuation)", flush=True)
print(f"G_PER_N={G_PER_N} | p_range=[50,20000]", flush=True)

for idx, gamma in ZEROS_71_100:
    print(f"\n{'='*70}", flush=True)
    print(f"γ{idx} (zero #{idx}) = {gamma:.6f}", flush=True)
    print(f"{'='*70}", flush=True)

    candidates = find_resonant_prime(gamma)
    if not candidates:
        print(f"  No candidates found!", flush=True)
        failed.append((idx, gamma, "no candidates"))
        continue

    print(f"Top-8 candidates in [50, 20000]:", flush=True)
    for err, p, j, geff in candidates:
        print(f"  p={p:>6}, j={j:>4}, γ_eff={geff:.6f}, err={err:.6f}", flush=True)

    winner = None
    winner_label = None
    for err, p, j, gamma_eff in candidates:
        N_list = N_SEQ_LARGE if p >= 5000 else N_SEQ_SMALL
        results = stability_test(p, N_list)
        label = classify_stability(results, p)
        row = f"  p={p:>6} j={j:>4} err={err:.6f}:"
        for N, log_ev0, delta in results:
            row += f"  N{N}={log_ev0:.4f}"
        row += f"  → {label}"
        print(row, flush=True)
        if winner is None and "STABLE" in label:
            winner = (err, p, j, gamma_eff, label)
            winner_label = label

    if winner is not None:
        err, p, j, gamma_eff, lbl = winner
        print(f"\nCONFIRMED γ{idx}: j={j}, lsq={p}, err={err:.6f}, {lbl}", flush=True)
        confirmed.append((idx, gamma, p, j, err, lbl))
    else:
        print(f"\nFAILED γ{idx}: no stable prime found in [50,20000]", flush=True)
        failed.append((idx, gamma, "no stable prime"))

print(f"\n{'='*70}", flush=True)
print(f"SUMMARY: {len(confirmed)}/30 confirmed, {len(failed)} failed", flush=True)
for idx, gamma, p, j, err, lbl in confirmed:
    print(f"  γ{idx:>3}: p={p}, j={j}, err={err:.6f}, {lbl}", flush=True)
if failed:
    print(f"FAILED:", flush=True)
    for item in failed:
        print(f"  γ{item[0]}: {item[2]}", flush=True)
