"""
Generic RMT-universality test suite — reusable across ANY ordered numeric
sequence (spectra, eigenvalues, gaps, arithmetic sequences).

Reuses the reference constants/curves already validated against zeta zeros
in ../../nn_ratio_exact.py, ../../exact_gue.py, ../../check_gue_variance.py,
../../pair_correlation_fixed.py, ../../analyze.py — this file generalizes
those one-off zeta-specific scripts into a candidate-agnostic pipeline.

Reference <r> constants (Atas et al. 2013, exact numerics, not the cruder
3x3-surmise approximation):
    Poisson : 0.38629  (= 2*ln(2) - 1)
    GOE (b1): 0.5307
    GUE (b2): 0.6027
    GSE (b4): 0.6769
"""

import math
import numpy as np
from scipy import stats

POISSON_R = 2 * math.log(2) - 1  # 0.386294
GOE_R = 0.5307
GUE_R = 0.6027
GSE_R = 0.6769


# =============================================================================
# 1. <r> ratio statistic — NO unfolding needed, run this first on everything
# =============================================================================

def mean_ratio(seq):
    """
    Oganesyan-Huse / Atas et al. ratio statistic.
    seq: sorted array-like of levels (raw, unfolded not required).
    Returns dict: mean, stderr, n, closest ensemble label.
    """
    x = np.asarray(sorted(seq), dtype=float)
    s = np.diff(x)
    s = s[s > 0]  # guard against duplicate levels
    if len(s) < 2:
        return {"mean": float("nan"), "stderr": float("nan"), "n": 0, "closest": None}
    r = np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])
    mean = float(np.mean(r))
    stderr = float(np.std(r, ddof=1) / math.sqrt(len(r)))

    refs = {"Poisson": POISSON_R, "GOE": GOE_R, "GUE": GUE_R, "GSE": GSE_R}
    closest = min(refs, key=lambda k: abs(refs[k] - mean))
    return {
        "mean": mean,
        "stderr": stderr,
        "n": len(r),
        "closest": closest,
        "sigma_from_gue": (mean - GUE_R) / stderr if stderr > 0 else float("inf"),
        "sigma_from_poisson": (mean - POISSON_R) / stderr if stderr > 0 else float("inf"),
    }


# =============================================================================
# 2. Unfolding — remove the smooth trend so mean spacing = 1
# =============================================================================

def unfold_poly(seq, degree=7):
    """
    Global polynomial unfolding: fit N(x) ~ polynomial(x) to the staircase
    counting function (x_i, i), then unfolded_i = fitted(x_i).
    Good when the sequence has a genuinely smooth global density trend
    (e.g. zeta zeros, most L-function zeros).
    """
    x = np.asarray(sorted(seq), dtype=float)
    n = len(x)
    idx = np.arange(1, n + 1, dtype=float)
    coeffs = np.polyfit(x, idx, degree)
    fitted = np.polyval(coeffs, x)
    # enforce monotonicity (polynomial fits can wiggle) by isotonic cleanup
    fitted = np.maximum.accumulate(fitted)
    unfolded = fitted
    spacings = np.diff(unfolded)
    return spacings[spacings > 0]


def unfold_local(seq, window=None):
    """
    Local-average unfolding: data-agnostic, no assumed functional form.
    Each raw spacing is divided by the local mean spacing over a symmetric
    window of neighboring spacings, then cumulatively summed so the result
    has mean spacing 1. Standard for sequences without a known analytic
    density (class numbers, graph eigenvalues, ideal norms, ...).
    """
    x = np.asarray(sorted(seq), dtype=float)
    raw = np.diff(x)
    raw = raw[raw > 0]
    n = len(raw)
    if window is None:
        window = max(15, n // 50)
    # local mean via centered moving average (edge-padded)
    kernel = np.ones(window) / window
    padded = np.pad(raw, (window // 2, window - window // 2 - 1), mode="edge")
    local_mean = np.convolve(padded, kernel, mode="valid")
    unfolded_spacings = raw / local_mean
    return unfolded_spacings


def unfold_rank(seq):
    """
    Trivial unfolding by construction: rescale so mean spacing is exactly 1
    globally (no local structure removed). Use only as a fallback / sanity
    check — this will NOT remove a nonuniform density trend, so P(s) and
    pair-correlation results from this method should be treated cautiously
    unless the sequence's density is already known-uniform.
    """
    x = np.asarray(sorted(seq), dtype=float)
    raw = np.diff(x)
    raw = raw[raw > 0]
    return raw / raw.mean()


# =============================================================================
# 3. Spacing distribution P(s)
# =============================================================================

def gue_wigner_pdf(s):
    return (32.0 / math.pi ** 2) * s ** 2 * np.exp(-4.0 * s ** 2 / math.pi)


def goe_wigner_pdf(s):
    return (math.pi / 2.0) * s * np.exp(-math.pi * s ** 2 / 4.0)


def poisson_pdf(s):
    return np.exp(-s)


def spacing_pdf(unfolded_spacings, bins=50, s_max=None):
    s = np.asarray(unfolded_spacings)
    if s_max is None:
        s_max = min(4.0, float(np.percentile(s, 99)))
    hist, edges = np.histogram(s, bins=bins, range=(0, s_max), density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, hist


def ks_against_ensembles(unfolded_spacings):
    """
    KS test of unfolded spacings against GUE / GOE / Poisson Wigner-surmise
    CDFs (closed-form for GOE/Poisson; GUE via numerical integration of the
    Wigner surmise — adequate for classification, not a high-precision claim).
    Returns dict of ensemble -> (D statistic, p-value).
    """
    s = np.asarray(unfolded_spacings)

    def gue_cdf(x):
        # CDF of (32/pi^2) s^2 exp(-4 s^2/pi) has a closed form via erf
        a = 4.0 / math.pi
        return stats.gamma.cdf(a * x ** 2, a=1.5) if x.ndim == 0 else np.array(
            [stats.gamma.cdf(a * xx ** 2, a=1.5) for xx in x]
        )

    def goe_cdf(x):
        return 1 - np.exp(-math.pi * x ** 2 / 4.0)

    def poisson_cdf(x):
        return 1 - np.exp(-x)

    out = {}
    for name, cdf in [("GUE", gue_cdf), ("GOE", goe_cdf), ("Poisson", poisson_cdf)]:
        d, p = stats.kstest(s, cdf)
        out[name] = {"D": float(d), "p": float(p)}
    return out


# =============================================================================
# 4. Pair correlation R2(x) vs GUE form 1 - (sin(pi x)/(pi x))^2
# =============================================================================

def pair_correlation(unfolded_cumsum_or_spacings, is_spacings=True, max_neighbors=20,
                      max_dist=5.0, bins=100):
    """
    unfolded_cumsum_or_spacings: either the unfolded spacings (is_spacings=True,
    default) or already-cumulative unfolded positions (is_spacings=False).
    """
    if is_spacings:
        positions = np.concatenate([[0.0], np.cumsum(unfolded_cumsum_or_spacings)])
    else:
        positions = np.asarray(unfolded_cumsum_or_spacings)
    N = len(positions)
    bin_width = max_dist / bins
    histogram = np.zeros(bins)

    for i in range(N):
        for j in range(i + 1, min(i + max_neighbors + 1, N)):
            diff = positions[j] - positions[i]
            if diff >= max_dist:
                break
            idx = int(diff / bin_width)
            if 0 <= idx < bins:
                histogram[idx] += 1

    centers = (np.arange(bins) + 0.5) * bin_width
    R2_emp = histogram / (N * bin_width)
    with np.errstate(divide="ignore", invalid="ignore"):
        sinc = np.sin(math.pi * centers) / (math.pi * centers)
    R2_gue = 1 - sinc ** 2
    return centers, R2_emp, R2_gue, N


# =============================================================================
# 5. Number variance / spectral rigidity Sigma^2(L)
# =============================================================================

def number_variance(unfolded_spacings, L_values=None, n_windows=2000, seed=0):
    """
    unfolded_spacings already has mean 1 (density-1 process).
    Sigma^2(L) for GUE ~ (2/pi^2)(ln(2 pi L) + gamma + 1 - pi^2/8)
    Sigma^2(L) for Poisson = L

    IMPORTANT: window start points must be drawn from the CONTINUUM,
    independent of the spectrum positions. Anchoring a window's left edge
    exactly on a spectrum point conditions on a point being there, which
    for a level-repulsive process (GUE) systematically suppresses the count
    right after the anchor and produces a spuriously flat/saturating
    Sigma^2(L) — caught 2026-07-18 when zeta-zero control showed Sigma^2(L)
    saturating at ~0.25 instead of growing like log(L).
    """
    positions = np.concatenate([[0.0], np.cumsum(unfolded_spacings)])
    N = len(positions)
    if L_values is None:
        L_values = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0]

    lo, hi = positions[N // 8], positions[7 * N // 8]
    rng = np.random.default_rng(seed)
    results = []
    euler_gamma = 0.5772156649015329
    for L in L_values:
        starts = rng.uniform(lo, hi - L, size=n_windows)
        # count of positions in [t0, t0+L) via searchsorted (positions sorted)
        lo_idx = np.searchsorted(positions, starts, side="left")
        hi_idx = np.searchsorted(positions, starts + L, side="left")
        counts = hi_idx - lo_idx
        var = float(np.var(counts))
        gue_var = (2 / math.pi ** 2) * (math.log(2 * math.pi * L) + euler_gamma + 1 - math.pi ** 2 / 8) if L > 0 else 0
        poisson_var = L
        results.append({"L": L, "empirical_variance": var, "gue_variance": gue_var,
                         "poisson_variance": poisson_var})
    return results


# =============================================================================
# 6. One-shot candidate runner
# =============================================================================

def run_full_suite(name, seq, unfold_method="local", n_min=200, notes=""):
    """
    Run the full pipeline on a candidate sequence and return a structured
    result dict. seq must have >= n_min distinct values or the deep tests
    are skipped (only <r> is reported).
    """
    x = np.asarray(sorted(set(np.asarray(seq, dtype=float).tolist())))
    result = {"name": name, "n": len(x), "notes": notes}

    r = mean_ratio(x)
    result["mean_ratio"] = r

    if len(x) < n_min:
        result["verdict"] = "insufficient N for deep tests (<r> only)"
        return result

    if unfold_method == "poly":
        unfolded = unfold_poly(x)
    elif unfold_method == "rank":
        unfolded = unfold_rank(x)
    else:
        unfolded = unfold_local(x)

    result["unfold_method"] = unfold_method
    result["unfolded_mean"] = float(np.mean(unfolded))
    result["unfolded_var"] = float(np.var(unfolded))

    centers, pdf = spacing_pdf(unfolded)
    result["spacing_pdf"] = {"centers": centers.tolist(), "pdf": pdf.tolist()}

    ks = ks_against_ensembles(unfolded)
    result["ks"] = ks

    pc_centers, pc_emp, pc_gue, N = pair_correlation(unfolded)
    skip = 5
    pc_mad = float(np.mean(np.abs(pc_emp[skip:] - pc_gue[skip:])))
    result["pair_correlation_mad"] = pc_mad

    nv = number_variance(unfolded)
    result["number_variance"] = nv

    # Decision rule: multiple independent signals must agree on GUE
    votes = []
    votes.append(r["closest"] == "GUE" and abs(r["sigma_from_gue"]) < 3)
    votes.append(min(ks, key=lambda k: ks[k]["D"]) == "GUE")
    votes.append(pc_mad < 0.08)
    nv5 = next((v for v in nv if v["L"] == 5.0), None)
    if nv5:
        votes.append(abs(nv5["empirical_variance"] - nv5["gue_variance"]) < 0.5 * nv5["poisson_variance"])
    n_votes = sum(bool(v) for v in votes)
    result["gue_votes"] = f"{n_votes}/{len(votes)}"
    result["verdict"] = "SHARES GUE PATTERN" if n_votes >= 3 else (
        "PARTIAL / AMBIGUOUS" if n_votes == 2 else "NOT GUE-like")

    return result


def summarize(result):
    r = result["mean_ratio"]
    line = f"{result['name']:40s} N={result['n']:>8d}  <r>={r['mean']:.4f}±{r['stderr']:.4f}  closest={r['closest']:6s}"
    if "verdict" in result:
        line += f"  verdict={result['verdict']}"
    return line
