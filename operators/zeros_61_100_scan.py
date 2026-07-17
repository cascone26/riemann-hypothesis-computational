#!/usr/bin/env python3
"""
zeros_61_100_scan.py

CvS spectral operator scan for Riemann zeros γ₆₁–γ₁₀₀.
Runs sequentially; designed to chain after zeros_51_60_scan.py completes.
Zero values from mpmath zetazero() at 25 decimal places.
"""
import math, numpy as np, scipy.linalg, sys
from scipy.special import digamma

TWO_PI = 2 * math.pi
LOGPI = math.log(math.pi)
G_PER_N = 900

ZEROS_61_100 = [
    (61, 165.537069187900),
    (62, 167.184439978175),
    (63, 169.094515415569),
    (64, 169.911976479412),
    (65, 173.411536519592),
    (66, 174.754191523366),
    (67, 176.441434297710),
    (68, 178.377407776100),
    (69, 179.916484020257),
    (70, 182.207078484366),
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


def prime_candidates(gamma, p_min=50, p_max=20000, top_n=8, err_threshold=None):
    cands = []
    for p in primes30000:
        if p < p_min or p > p_max:
            continue
        log_p = math.log(p)
        j = round(gamma * log_p / TWO_PI)
        if j < 1:
            continue
        ge = TWO_PI * j / log_p
        err = abs(ge - gamma)
        if err_threshold is not None and err >= err_threshold:
            continue
        cands.append((err, j, p, ge))
    cands.sort()
    return cands[:top_n]


def test_candidate(p, gamma, nvals=None):
    if nvals is None:
        n_max = 320 if p < 5000 else 256
        nvals = [64, 128, 192, n_max]
    row = {}
    for N in nvals:
        e = compute_ev0(p, N)
        row[N] = math.log10(abs(e)) if abs(e) > 0 else -99
    n_last = nvals[-1]
    delta = row[n_last] - row[64]
    stable = (row[n_last] < -10) and (delta < 2.5)
    all_stable = all(row[N] < -10 for N in nvals)
    vals = "  ".join(f"N{N}={row[N]:>7.4f}" for N in nvals)
    if all_stable:
        verdict = "STABLE***"
    elif stable:
        verdict = f"STABLE(1side,Δ={delta:+.2f})"
    else:
        verdict = f"div(Δ={delta:+.2f})"
    return row, delta, stable, all_stable, verdict, vals


def scan_zero(k, gamma):
    print(f"\n{'='*70}")
    print(f"γ{k} (zero #{k}) = {gamma:.6f}")
    print(f"{'='*70}")
    sys.stdout.flush()

    cands = prime_candidates(gamma, p_min=50, p_max=20000, top_n=8)
    print(f"Top-8 candidates in [50, 20000]:")
    for err, j, p, ge in cands:
        print(f"  p={p:6d}, j={j:4d}, γ_eff={ge:.6f}, err={err:.6f}")

    confirmed = None
    for err, j, p, ge in cands:
        row, delta, stable, all_stable, verdict, vals = test_candidate(p, gamma)
        print(f"  p={p:6d} j={j:4d} err={err:.6f}: {vals}  → {verdict}", flush=True)
        if stable and confirmed is None:
            confirmed = (j, p, err, ge, verdict)

    if confirmed is not None:
        j, p, err, ge, verdict = confirmed
        print(f"\nCONFIRMED γ{k}: j={j}, lsq={p}, err={err:.6f}, {verdict}")
        return confirmed

    print(f"\nRound 1 failed. Extending error threshold to err<0.01 in [50, 20000]...")
    cands_wide = prime_candidates(gamma, p_min=50, p_max=20000, top_n=20, err_threshold=0.01)
    already_tested = {p for _, _, p, _ in cands}
    new_cands = [(err, j, p, ge) for err, j, p, ge in cands_wide if p not in already_tested]
    if new_cands:
        print(f"  Testing {len(new_cands)} new candidates with 0.001 ≤ err < 0.01:")
        for err, j, p, ge in new_cands[:12]:
            row, delta, stable, all_stable, verdict, vals = test_candidate(p, gamma)
            print(f"  p={p:6d} j={j:4d} err={err:.6f}: {vals}  → {verdict}", flush=True)
            if stable and confirmed is None:
                confirmed = (j, p, err, ge, verdict)
                break

    if confirmed is not None:
        j, p, err, ge, verdict = confirmed
        print(f"\nCONFIRMED γ{k}: j={j}, lsq={p}, err={err:.6f}, {verdict}")
        return confirmed

    print(f"\nRound 2 failed. Extending prime range to [20000, 30000]...")
    cands_large = prime_candidates(gamma, p_min=20000, p_max=30000, top_n=8)
    if cands_large:
        print(f"  Top-8 candidates in [20000, 30000]:")
        for err, j, p, ge in cands_large:
            row, delta, stable, all_stable, verdict, vals = test_candidate(p, gamma)
            print(f"  p={p:6d} j={j:4d} err={err:.6f}: {vals}  → {verdict}", flush=True)
            if stable and confirmed is None:
                confirmed = (j, p, err, ge, verdict)
                break

    if confirmed is not None:
        j, p, err, ge, verdict = confirmed
        print(f"\nCONFIRMED γ{k}: j={j}, lsq={p}, err={err:.6f}, {verdict}")
        return confirmed

    print(f"\nγ{k} NOT CONFIRMED after 3 rounds. Needs further investigation.")
    return None


print("CvS Spectral Operator — Zeros 61-100 Scan")
print(f"Stability criterion: N_max < -10 AND delta(N_max−N64) < +2.5")
print(f"G_PER_N={G_PER_N}, prime sieve up to 30000")
print()

results = []
for k, gamma in ZEROS_61_100:
    result = scan_zero(k, gamma)
    results.append((k, gamma, result))
    sys.stdout.flush()

print("\n\n" + "=" * 70)
print("SUMMARY — Zeros 61-100")
print("=" * 70)
confirmed_count = 0
for k, gamma, result in results:
    if result is not None:
        j, p, err, ge, verdict = result
        print(f"  γ{k:3d} = {gamma:.6f}: CONFIRMED j={j}, lsq={p}, err={err:.6f}, {verdict}")
        confirmed_count += 1
    else:
        print(f"  γ{k:3d} = {gamma:.6f}: NOT CONFIRMED")

print(f"\nConfirmed: {confirmed_count}/40 (running total: {60+confirmed_count}/100)")
print("\nDone.")
