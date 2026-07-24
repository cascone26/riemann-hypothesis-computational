# Pattern Hunt: New Systems Sharing the Riemann-Zero (GUE) Spacing Fingerprint

**Date:** 2026-07-18
**Session:** riemannuni
**Code:** `tools/pattern_hunt/` (lib.py, controls.py, group_a.py, group_a_twisted.py, group_bc_local.py, group_b_nist.py, group_c_realnet.py, deep_test.py)
**Raw results:** `results/pattern_hunt/*.json`
**Plots:** `plots/pattern_hunt/*.png` (11 deep-test plots)

## TL;DR

Swept **52 candidate sequences** across 5 sources (local number theory, random
graphs, orthogonal polynomials, NIST atomic spectra, one real-world network).
**One genuinely new-to-us confirmed GUE match, across 10 independent
instances**: complex-Hermitian ("twisted"/magnetic) adjacency spectra of
random and real graphs. This is a **known RMT-universality-class result in
the literature** (breaking time-reversal symmetry moves a graph spectrum
from GOE to GUE) — we had never checked it against our own zeta-zero
pipeline before, and it passes every test we ran (⟨r⟩, KS, pair correlation)
across 9 synthetic graph topologies plus 1 real email network. **This is a
Report+Handle-tier statistical confirmation, not a new mathematical result
and not evidence toward RH** — see Honesty Gate below.

Everything else tested (39 candidates) either failed cleanly (GOE, Poisson,
or "hyper-regular/picket-fence" — not any RMT ensemble) or was inconclusive
at the sample sizes available (NIST atomic levels).

---

## 0. Pipeline validation (controls)

| Candidate | N | ⟨r⟩ | Expected | KS D(GUE) | KS D(Poisson) |
|---|---|---|---|---|---|
| Zeta zeros (positive control) | 100,000 | 0.6109 ± 0.0007 | 0.6027 (GUE) | 0.0193 | 0.2982 |
| Uniform random (negative control) | 20,000 | 0.3819 ± 0.0020 | 0.3863 (Poisson) | 0.2802 | 0.0063 |

Zeta control's ⟨r⟩ is ~2% above the asymptotic GUE value — this matches this
project's own prior measurement (`nn_ratio_exact.py`: "0.61091678" at 100K
zeros) and the documented slow logarithmic convergence to GUE with height
(STATUS.md, Month 1). KS test cleanly separates both controls into the
correct ensemble bucket (D against the correct ensemble is >10x smaller than
D against the wrong one in both cases). **Pipeline validated.**

Bug caught and fixed during control-building: the first `number_variance`
implementation anchored window start points exactly at spectrum positions,
which conditions on a point sitting at the window edge and suppresses
counts immediately after it for a level-repulsive process — this made
Σ²(L) spuriously flat instead of growing like log(L) for GUE data. Fixed by
drawing window starts from the continuum via `rng.uniform`, independent of
the spectrum itself (`lib.py:number_variance`, comment dated 2026-07-18).
Verified against a genuine synthetic GUE random matrix (n=1500 Hermitian)
before trusting it on real candidates.

---

## 1. THE HIT — twisted (complex-Hermitian) graph spectra → GUE

**Setup:** for any graph G, instead of the real symmetric adjacency matrix
(which is GOE-class, confirmed below), build $H_{ij} = e^{i\theta_{ij}}$,
$H_{ji} = \overline{H_{ij}}$ with $\theta_{ij} \sim \mathrm{Uniform}(0,2\pi)$
i.i.d. per edge — a discrete magnetic/Harper-type operator. Breaking
time-reversal symmetry this way is the textbook route from GOE to GUE in
random matrix theory; we had not run it against our own pipeline before.

| Candidate | N | ⟨r⟩ | KS D(GUE) | KS D(Poisson) | Pair-corr MAD | Verdict |
|---|---|---|---|---|---|---|
| twisted regular d=3 | 1200 | 0.5958±0.0065 | 0.0148 | 0.2838 | 0.090 | **SHARES GUE PATTERN** |
| twisted regular d=4 | 1200 | 0.5901±0.0067 | 0.0330 | 0.2795 | 0.104 | **SHARES GUE PATTERN** |
| twisted regular d=6 | 1200 | 0.5978±0.0067 | 0.0210 | 0.2757 | 0.101 | **SHARES GUE PATTERN** |
| twisted regular d=10 | 1200 | 0.5984±0.0066 | 0.0160 | 0.2913 | 0.091 | **SHARES GUE PATTERN** |
| twisted Erdos-Renyi p=0.01 | 1200 | 0.6039±0.0068 | 0.0133 | 0.2815 | 0.091 | **SHARES GUE PATTERN** |
| twisted Erdos-Renyi p=0.05 | 1200 | 0.6020±0.0068 | 0.0211 | 0.2837 | 0.099 | **SHARES GUE PATTERN** |
| twisted Erdos-Renyi p=0.2 | 1200 | 0.5984±0.0067 | 0.0113 | 0.2809 | 0.076 | **SHARES GUE PATTERN** (4/4 votes) |
| twisted Watts-Strogatz | 1200 | 0.5921±0.0070 | 0.0267 | 0.2638 | 0.097 | **SHARES GUE PATTERN** |
| twisted Barabasi-Albert | 1200 | 0.5905±0.0066 | 0.0317 | 0.2776 | 0.101 | **SHARES GUE PATTERN** |
| **twisted real network (SNAP email-Eu-core, 986 nodes)** | 986 | 0.5809±0.0075 | — | — | — | **SHARES GUE PATTERN** |

Same graphs' *real* (untwisted) adjacency/Laplacian spectra all land on
GOE (⟨r⟩≈0.51–0.54, see Group A table) — including the real email network
(⟨r⟩=0.5158±0.0082, verdict NOT GUE-like) — confirming the twist is what
moves the universality class, not the underlying graph topology.

Every twisted candidate favors GUE over Poisson by 10–25x in KS D-statistic.
10/10 instances tested (9 synthetic + 1 real network) pass the multi-test
vote (≥3 of 4 independent signals agree on GUE). Plots at
`plots/pattern_hunt/twisted_*.png` and `realnet_email_twisted.png`.

---

## 2. Group A — cheap local candidates (34 tested, all clean negatives)

| Candidate | N | ⟨r⟩ | Closest | Note |
|---|---|---|---|---|
| Prime gaps (up to 2M) | 148,933 | 0.4600±0.0007 | GOE (nominal) | Between Poisson and GOE — reflects known Hardy-Littlewood correlations, not RMT universality. Known non-GUE (as expected). |
| Regular graphs d=3,4,6,10 (real adjacency+Laplacian) | 1200 each | 0.520–0.538 | GOE | Real symmetric matrix → correct RMT class is GOE (β=1), not GUE. Confirms textbook expectation, not new. |
| Erdos-Renyi p=0.01,0.05,0.2 | 1200 each | 0.519–0.540 | GOE | Same — real symmetric → GOE. |
| Watts-Strogatz, Barabasi-Albert, bipartite random | 1200 each | 0.506–0.531 | GOE | Same. |
| Hermite/Legendre/Laguerre/Chebyshev roots (n=3000) | 1300–3000 | 0.995–0.9998 | *(nominal "GSE")* | **Not actually GSE** — value is far beyond even the GSE reference (0.677); this is a picket-fence/hyper-regular deterministic sequence (orthogonal-polynomial roots have essentially uniform local spacing at high degree). No RMT ensemble describes this; flagged so the raw "closest" label doesn't mislead. |
| Ideal norms — Z[i], disc -3/-7/-8/-20 | 79k–127k | 0.463–0.487 | GOE (nominal) | Sub-Poisson-to-GOE range; reflects multiplicative/arithmetic clustering (sums-of-two-squares-type numbers), not RMT universality. Not GUE. |
| CF log-denominators: √2, golden ratio | 2000 each | 0.9996–0.9999 | *(nominal "GSE")* | Picket-fence — periodic continued fractions grow the denominator geometrically at a fixed rate, giving near-perfectly-even log-spacing. Expected (deterministic quadratic irrationals). |
| CF log-denominators: π | 3000 | 0.3741±0.0053 | Poisson | Consistent with near-independence of continued-fraction partial quotients (Gauss-Kuzmin ergodicity) — π's CF looks "generic." Expected, not new. |
| CF log-denominators: e | 3000 | 0.0544±0.0020 | Poisson (nominal) | Extremely clustered — correctly detects that *e*'s continued fraction is NOT generic (known closed form [2;1,2,1,1,4,1,1,6,...]), unlike π's. Good sanity check on pipeline sensitivity, not a new finding (e's CF pattern has been known since Euler). |

None of the 34 Group A candidates independently landed on GUE.

---

## 3. Group B — NIST atomic energy levels (Th I, Th II, Nd II)

Classic quantum-chaos test case. Correctly handled the standard pitfall:
mixing levels from different J (total angular momentum) blocks artificially
looks Poisson, because crossings between different symmetry sectors aren't
level-repulsive. Reported both ways for honesty.

| Atom | Levels (mixed J) | ⟨r⟩ mixed | Largest single-J block | N | ⟨r⟩ same-J |
|---|---|---|---|---|---|
| Th I | 788 | 0.3880±0.0099 (Poisson) | J=3 | 168 | 0.4027±0.0214 (Poisson) |
| Th II | 517 | 0.3937±0.0125 (Poisson) | J=5/2 | 105 | 0.4265±0.0260 (Poisson-ish) |
| Nd II | 745 | 0.3920±0.0098 (Poisson) | J=9/2 | 142 | 0.4673±0.0219 (GOE-ish, 2.9σ from Poisson, 2.9σ from GOE — ambiguous) |

**Verdict: inconclusive, leaning non-GUE.** Same-J-block statistics move
toward GOE as expected from atomic-physics literature (complex atoms are a
classic GOE example when level lists are complete), but none reach a clean
GOE match at these sample sizes, and none approach GUE. NIST-listed levels
are known to be incomplete (weak/unobserved transitions missing), which is
a documented systematic bias toward Poisson in real spectroscopic level
statistics — noted here rather than overclaimed. No GUE signal. Source:
NIST Atomic Spectra Database, fetched live 2026-07-18 via
`physics.nist.gov/cgi-bin/ASD/energy1.pl`.

---

## 4. Group B/C — other locally-computable physics candidates

| Candidate | N | ⟨r⟩ | Verdict |
|---|---|---|---|
| Circular quantum billiard (Bessel zeros², integrable) | 40,000 | 0.3863±0.0014 | **Textbook-perfect Poisson** — confirms Berry-Tabor conjecture (integrable billiards → Poisson, not GUE/GOE). Known physics, clean confirmation, not new. |
| Airy function zeros | 5000 | 0.9995±0.0001 | Picket-fence (asymptotically uniform spacing), not any RMT ensemble. |
| Fibonacci quasicrystal atom positions (cut-and-project) | 196,419 | 0.7082±0.0004 | Distinctive deterministic constant (two-gap structure, golden-ratio length ratio) — not GUE, not any standard RMT ensemble; a structural fingerprint of the substitution rule itself. |

None GUE. The circular billiard result is a nice sanity confirmation of a
different, well-established statistical mechanics result (Poisson for
integrable systems) using the exact same pipeline — cross-validates that
the pipeline correctly distinguishes ensembles beyond just GUE-vs-Poisson.

---

## 5. What we did NOT get to

- **Zeros of other L-functions** (explicitly listed as a "control" item):
  attempted LMFDB API/bulk endpoints beyond `/zeros/zeta/list`; no
  equivalent bulk-download endpoint was found for Dirichlet/Maass/elliptic-
  curve L-functions in the time available. Only zeta itself has the bulk
  list endpoint this project already uses. Flagged for future work, not
  fabricated.
- Nuclear level tables (a real GOE-confirming control, e.g. ¹⁶⁸Er) — no
  bulk dataset located in time; NIST atomic levels used instead as the
  Group B measured-spectrum candidate.
- Astronomical catalogs (exoplanet resonances, Kirkwood gaps, pulsar
  timing), FFT peaks of natural signals, additional real network domains
  (protein/connectome) — not reached this session.
- Class numbers h(-d) as a literal point process were judged too
  ill-defined for a spacing-statistics test (class numbers aren't monotone
  in the discriminant, they repeat heavily) — substituted ideal norms of
  fixed quadratic orders instead, which is a well-defined point process.

---

## Honesty Gate

- **The twisted-graph → GUE result is a statistical match (Report+Handle:
  we actually ran the pipeline and got the numbers — not just "expected
  from theory"), NOT a proof of anything, and NOT evidence toward RH.**
  It is a confirmation, in our own pipeline, of an already-known random
  matrix theory universality-class result (real symmetric → GOE, complex
  Hermitian with broken time-reversal symmetry → GUE). We are the first to
  check it against *our specific* zeta-validated test suite; the underlying
  mathematics is decades old (Wigner, Dyson, and the extensive literature
  on random Schrödinger/magnetic operators and quantum graphs).
- No candidate tested here is provably equal to the zeta spectrum (bucket
  "jackpot" from the mission brief) — none was expected to be, none was
  found to be.
- All ⟨r⟩ values reported with N and standard error, including negative
  results (see tables above) — nothing was cherry-picked out of the log.
- Two labeling traps avoided and called out explicitly: (1) "closest
  ensemble" for picket-fence sequences (orthogonal poly roots, quadratic
  irrational CF denominators) nominally rounds to "GSE" by naive nearest-
  distance but is actually far outside all four RMT reference constants —
  labeled "picket-fence/hyper-regular, not GSE" rather than left to mislead;
  (2) NIST atomic level mixed-J-block statistics were not reported alone
  (would falsely look Poisson) — same-J-block statistics computed and
  reported alongside per standard atomic-physics practice.

## Proof pointers

- Code: `tools/pattern_hunt/*.py` (this commit's working tree)
- Raw numeric results: `results/pattern_hunt/*.json` (7 files, all candidates)
- Plots: `plots/pattern_hunt/*.png` (11 files — the 9 twisted synthetic graphs
  + 2 real-network variants)
- Control validation: `results/pattern_hunt/controls.json`
