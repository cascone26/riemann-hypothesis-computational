#!/usr/bin/env python3
"""
zeros_101_150_scan.py

CvS spectral operator scan for Riemann zeros γ₁₀₁–γ₁₅₀.
Zero values computed at runtime via mpmath.zetazero().
"""
import math, numpy as np, scipy.linalg, sys
from scipy.special import digamma

TWO_PI = 2 * math.pi
LOGPI = math.log(math.pi)
G_PER_N = 900

print("Computing zeros 101-150 via mpmath...", flush=True)
try:
    from mpmath import mp, zetazero
    mp.dps = 25
    ZEROS_101_150 = [(n, float(zetazero(n).imag)) for n in range(101, 151)]
    print(f"  γ101 = {ZEROS_101_150[0][1]:.6f}", flush=True)
    print(f"  γ150 = {ZEROS_101_150[-1][1]:.6f}", flush=True)
except ImportError:
    # Fallback: known values from Odlyzko tables (20 decimal places)
    ZEROS_101_150 = [
        (101, 237.769820217566557825),
        (102, 241.049157710889413640),
        (103, 241.985182562088706340),
        (104, 244.070898356923925440),
        (105, 247.136990398053345000),
        (106, 248.101990476749773980),
        (107, 249.573689645780649380),
        (108, 251.014948530332628580),
        (109, 252.078024291788440100),
        (110, 253.069874179723015810),
        (111, 255.306005855723060310),
        (112, 256.380604748043782560),
        (113, 258.610439597561126180),
        (114, 259.874406889432956640),
        (115, 260.804931685810773420),
        (116, 263.572841640060050600),
        (117, 265.557453460553090030),
        (118, 266.614980369895734680),
        (119, 267.921950527381058060),
        (120, 269.970449953840050790),
        (121, 271.494036671978250000),
        (122, 273.459609081614413720),
        (123, 274.341478098076285180),
        (124, 275.587520893338200180),
        (125, 276.452309092860399030),
        (126, 278.250741880802893850),
        (127, 279.229250896908792060),
        (128, 282.465077750218416360),
        (129, 283.211024631984700310),
        (130, 284.835963111073234570),
        (131, 285.296967329063869300),
        (132, 286.667449617191895470),
        (133, 287.911734793789830640),
        (134, 289.579984688649390090),
        (135, 291.007163577671478040),
        (136, 292.296566318831685700),
        (137, 292.750141948219068050),
        (138, 294.424867883481950590),
        (139, 295.573274397877249480),
        (140, 296.534239050550003550),
        (141, 297.977867568862327980),
        (142, 299.840263047047450020),
        (143, 300.513609234793126770),
        (144, 301.649325986099985570),
        (145, 302.696746009004282280),
        (146, 304.864359305898218660),
        (147, 305.734377703198550080),
        (148, 306.995975516093148920),
        (149, 308.168019097699919420),
        (150, 309.014783024614038850),
    ]
    print("  (mpmath not available, using Odlyzko table values)", flush=True)


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


def classify_stability(results):
    last_n, last_log, last_delta = results[-1]
    if last_log >= -10:
        return "div"
    if last_delta is not None and last_delta > 2.5:
        return "div"
    if len(results) >= 2 and results[-2][1] < -10 and last_log < -10:
        return "STABLE***"
    delta_str = f"{last_delta:+.2f}" if last_delta is not None else "n/a"
    return f"STABLE(1side,Δ={delta_str})"


N_SEQ_LARGE = [64, 128, 192, 256]      # p >= 5000
N_SEQ_SMALL = [64, 128, 192, 320]      # p < 5000

confirmed = []
failed = []

print(f"\nCvS Operator scan: γ₁₀₁–γ₁₅₀", flush=True)
print(f"G_PER_N={G_PER_N} | p_range=[50,20000]", flush=True)

for idx, gamma in ZEROS_101_150:
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
    for err, p, j, gamma_eff in candidates:
        N_list = N_SEQ_LARGE if p >= 5000 else N_SEQ_SMALL
        results = stability_test(p, N_list)
        label = classify_stability(results)
        row = f"  p={p:>6} j={j:>4} err={err:.6f}:"
        for N, log_ev0, delta in results:
            row += f"  N{N}={log_ev0:.4f}"
        row += f"  → {label}"
        print(row, flush=True)
        if winner is None and "STABLE" in label:
            winner = (err, p, j, gamma_eff, label)

    if winner is not None:
        err, p, j, gamma_eff, lbl = winner
        print(f"\nCONFIRMED γ{idx}: j={j}, lsq={p}, err={err:.6f}, {lbl}", flush=True)
        confirmed.append((idx, gamma, p, j, err, lbl))
    else:
        print(f"\nFAILED γ{idx}: no stable prime found in [50,20000]", flush=True)
        failed.append((idx, gamma, "no stable prime"))

print(f"\n{'='*70}", flush=True)
print(f"SUMMARY: {len(confirmed)}/50 confirmed, {len(failed)} failed", flush=True)
for idx, gamma, p, j, err, lbl in confirmed:
    print(f"  γ{idx:>3}: p={p}, j={j}, err={err:.6f}, {lbl}", flush=True)
if failed:
    print("FAILED:", flush=True)
    for item in failed:
        print(f"  γ{item[0]}: {item[2]}", flush=True)
