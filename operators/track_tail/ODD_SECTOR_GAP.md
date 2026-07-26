# ODD SECTOR GAP — piece 3 cert is vacuous; sector never certified; machinery transfers
*(Hawk, 2026-07-24 roundtable. Receipts in this doc; measurements re-runnable via the scripts cited.)*

## 1. The finding

THEOREM_ASSEMBLY piece 3 ("Odd sector cert — DONE/Handle") rests on `arb_certify.py --parity odd`,
which certifies `modes = [1]`. But per LEDGER §192 (parity decoupling, Handle, machine-zero off-block
4.8e-39): **EVEN-v sector = odd-k sines, ODD-v sector = even-k sines.** k=1 is an odd sine index —
EVEN parity. The tool has no parity branch in its entries (`dterms`/`Q` identical for both flags;
only the mode list changes), so `--parity odd` certifies a 1×1 principal minor of the EVEN block.

**Numerical receipt:** the odd cert's pivot [0.0002395 ± 5.56e-8] equals A·Q(1,1) of the even-sector
engine exactly: A·Qraw(1,1) = 0.549 · 4.36325e-4 = **2.39543e-4** (float NL=20000; arb_certify's Q
omits the /A normalization that k0_hp's Qraw includes). Not a coincidence, an identity.

**Consequence:** the odd sector (even k) had ZERO Arb-certified entries before tonight, and the
"minimal viable theorem is assembly-ready" line in THEOREM_ASSEMBLY was wrong on this count
(independently of the known-open piece-2 tail). §192's odd-sector numbers (full λmin 9.644e-5 at
a=0.549) are float/Report-tier — the "generous margin, easy" comment in the cert is a heuristic
from a dead session, the exact failure mode this project keeps catching.

## 2. First-ever odd-sector (even-k) measurements — and full transfer of tonight's machinery

Float kernel, NL=20000 (`uc_closed_form_check.py` functions, k even):

- **Diagonals:** Q(2,2)=0.02162 (small corner!), Q(4,4)=0.4273, Q(6,6)=0.6709, Q(8,8)=1.5939,
  Q(26,26)=2.3435, Q(40,40)=3.1499, Q(50,50)=3.0492. Converges to the SAME density law
  log(kb/2π) − 0.3614·cos(kb·log2): at k=40/50 agreement to 0.7%/0.13%. Deviations at small k are
  the parity-dependent 1/k terms (cos(kπ)=+1 for even k flips boundary signs).
- **NOT diagonal-dominant at the corner:** |Q(2,4)|=0.0936 > Q(2,2)=0.0216; |Q(4,6)|=0.3623 ~ 85%
  of Q(4,4). The sector needs a genuine low-block LDL cert (like the even sector's), not Gershgorin.
- **Loewner form holds on even modes at machine precision:** Q[k,l] = (2kl/(A(k²−l²)))(ĥ(l)−ĥ(k))
  verified at rel. err 1.6e-13 .. 1.7e-15 on (2,4),(2,6),(4,6),(2,26),(16,40). SAME ĥ sequence.
- **u-law transfers with the parity sign flip:** for even m, Si(mπ) = π/2 − 1/(mπ), so
  u(m) = −1 − (4w_p/π)sin(θm) + 2/(π²m) − 0.162/m; measured resid·m = −0.162 constant to 3 digits
  (m = 26..250). u_c = −1 and β = 4w_p/π carry over EXACTLY; only the O(1/m) constants are
  parity-dependent (odd m: +1.517 net; even m: −0.162−0.203 net — to be pinned in the writeup).

So the odd sector inherits, with no new analysis beyond parity-sign bookkeeping: the Cauchy π/2
Hilbert bound, the β-channel, the rank-1 c/(kl) remainder, the diagonal density law, and the 4′
cut architecture. The theorem becomes: [even: certified low block + 4′] ⊕ [odd: certified low
block + mirrored 4′], composed by exact parity decoupling (§192, machine zero).

## 3. Status / actions

- `arb_certify_oddparity_k0.py` created (variant of arb_certify_k0.py; modes {2,4,...,K0}), fired at
  K0=24, NL=50000 → `results/arb_certify_oddparity_K24_nl50000.log`. **RESULT: CERTIFIED PD, exit 0 —
  all 12 pivots rigorously > 0** (mode 24: [1.48133 ± 3.78e-6] down to the corner mode 2:
  [1.12e-5 ± 3.96e-8], radius 280× inside the value). First rigorous odd-sector certification in the
  project. The k=2 corner survived LDL against the large neighbors exactly as the large-diagonal-first
  ordering intended. Piece 3 now has a real low-block artifact; the odd TAIL (K0_odd composition +
  mirrored 4′) remains open as scoped above.
- THEOREM_ASSEMBLY piece 3 row: correction blockquote added (do not delete the old row — provenance).
- **UPDATE (same session):** the odd composition table is DONE — `k0_hp_oddparity_nl24000.py`,
  `results/k0_hp_oddparity_nl24000_out.txt`, `results/Qodd52_a0549_nl24000_float.npy`. Composition
  HOLDS at the certified K0=24: λmin 1.888e-5 (whitened) vs near-null residual 1.312e-6 → headroom
  14.4× (minimal K0_odd=10 at 1.13× — too thin, avoid). Direct tail-window λmin 2.029 confirms the
  d₀ law. NORMALIZATION NOTE: Arb certs are unwhitened (×A vs these tables) — arb mode-2 pivot
  1.12e-5 ↔ whitened 2.04e-5 ≥ λmin 1.888e-5, interlacing-consistent.
- **§222 REPAIR SWEEP (fixed multiprime variant, modes {2,4,6,8}, NLERCH=50000,
  `results/multiprime_oddparity_sweep.log`):** TRUE odd sector CERTIFIED PD at a ∈ {0.55, 0.58,
  0.62, 0.64, 0.70}; at a ∈ {0.75, 0.80} NOT certified — mode-2 pivot interval contains 0
  ([±1.62e-7], [±6.82e-8]) with a monotone collapse trend (3.08e-5 → 2e-7 over 0.55→0.70), i.e.
  undetermined-at-this-precision approaching the p=5 criticality (a=½log5≈0.805) — the odd twin of
  §225's even-margin dip above 0.64. NET: §225's "both-parity reach [0.549, 0.64]" is RESTORED with
  a genuine odd sector behind it (odd actually reaches 0.70; even remains the binding constraint).
- **GATE BUG (found via this sweep):** `arb_certify_multiprime.py` prints "NOT certified" but exits 0
  — no failure exit code. Fixed in the `_oddparity` variant (sys.exit(1) on failure); the ORIGINAL
  still has the bug — flagged for the pieces-5/6 audit lane, same §207 failure class.
- Ledger entries due: §235 (census), §236 (K0=23 cert), §237 (Loewner + u_c=−1 + diagonal law),
  §238 (4′ architecture + no-split retirement), §239 (this gap + odd-sector transfer + cert result).
