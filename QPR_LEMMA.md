# The ‖Q_pr‖ Lemma — operator bound for the cross term of the p⊕r assembly (item 2; drafts §803)

**Goal.** For p in the certified probe span and r the class residual, bound |Q(p,r)| ≤ B·‖p‖‖r‖
(L² norms) with B an explicit certified constant. The assembly inequality (§801) then reads
λ_min‖p‖² − 2B‖p‖‖r‖ − ε_R‖r‖² ≥ 0, which holds for ‖r‖/‖p‖ below a computable threshold —
and the class definition (top-eigenspace of the double-limiting operator) supplies that smallness
via ε_R.

## Decomposition: B ≤ B_arch + B_pole + B_prime, each certified separately

All three via the same pattern: Q's sesquilinear pieces against L²-normalized band-limited
functions, using only (i) band-limit 4W (time support), (ii) the window's quartic Fejér envelope,
(iii) closed forms already in the apparatus.

**B_prime (easy, closed form).** |Σ_n w_n ĝ_{pr}(log n)| with w_n = Λ(n)/√n, and
|ĝ_{pr}(u)| = |(1/2π)∫ p̂(t) r̂(t) cos(ut) dt| ≤ (1/2π)‖p̂‖₂‖r̂‖₂ = (1/2π)‖p‖‖r‖ (Plancherel
convention as in the apparatus; factor 2 from the implemented 2Σ noted at assembly). So
B_prime = (1/π)·Σ_{log n<4W} Λ(n)/√n. At 4W=1.32: Σ = log2/√2 + log3/√3 = 1.1244 ⇒
**B_prime ≤ 0.3579** (exact arithmetic; per-width formula, grows ~e^{2W}).

**B_pole (rank-one, closed form).** pole(p,r) = 2·Re[p̂(i/2) r̂(i/2)] (apparatus convention).
For band-limited ψ with time support [−2W,2W]: ψ̂(i/2) = ∫ψ_time(u)e^{u/2}du ≤ e^{W}·‖ψ_time‖₁
≤ e^{W}·√(4W)·‖ψ‖₂ (Cauchy–Schwarz on the support). So B_pole ≤ 2·e^{2W}·4W·(1/2π)-normalized —
CONVENTION-SENSITIVE: derive against the exact pole_entry normalization in arb_fejer_cert2 before
certifying; the shape is B_pole = C(W)·e^{2W} with C explicit. At 4W=1.32 expect O(1-3).

**B_arch (the one careful piece).** arch(p,r) = (1/2π)∫ p̂ r̂ Ω(t) dt, Ω(t) = Re ψ₀(¼+it/2) − ln π
~ ln|t|. Split the line at the window slack [T₋, T₊] = [0, 2·C_max]:
  - inside: |∫_{|t|≤T₊}| ≤ (1/2π)·max_{|t|≤T₊}|Ω| · ‖p‖‖r‖ ≤ (ln T₊ + O(1))/2π ·‖p‖‖r‖ ≈
    (ln 88 + 2)/2π ≈ **1.03**·‖p‖‖r‖ at C_max=44 (digamma asymptotic bound as in
    arb_zerofree_cert3's header, certified there);
  - outside: both p̂ and r̂ carry the window envelope. p is in the probe span: |p̂(t)| ≤
    ‖coef‖₁·2/(W(t−44))² (§791's envelope, certified). r is in the class: same quartic envelope
    holds for the double-limited class by definition of the window limiting (STATE precisely in
    the class definition — this is where the class choice does work). Then
    |∫_{t>T₊}| ≤ (2/π)·(lnT₊+2)/(W²(T₊−44)³)-type tail — same integral family §791 certified;
    at T₊=88 it is ~1e-4, negligible.
  So **B_arch ≈ 1.03 + tails** at the anchor — dominated by the in-window digamma sup.

**Total at 4W=1.32, C_max=44: B ≲ 1.4 + B_pole(TBD convention).** Against λ_min = 0.1289 (§802)
the assembly needs ‖r‖/‖p‖ ≲ 0.1289/(2B) ≈ 0.04-0.05 — i.e. ε_R ≲ 2e-3 in eigenvalue terms.
The prolate plunge at capacity+2 is orders below that (measured §794 spectra: 4th eigenvalue
already ≥ 6e-2 vs kernel-tier 1e-6 — the gap is wide). The assembly CLOSES comfortably if ε_R is
certified anywhere below 1e-3.

## Status & next
- B_prime: exact, done above (certify in one line of arb).
- B_arch: derivation above; certify sup|Ω| on [0,88] (arb, one integral-free bound) + the §791-
  family tail. One short script.
- B_pole: pin the pole_entry convention first (read arb_fejer_cert2.pole_entry), then closed form.
- ε_R: route A = KRD constants (PDF-read pending, §798); route B (BETTER for the anchor) =
  per-instance: define the class as top-k eigenspace of the double-limiting operator, certify
  λ_{k+1} ≤ ε directly with the interval pencil on a padded probe set — our own §800 "S NOT PD at
  n=16" failure is this spectrum's plunge showing in the instrument; measure it as a curve instead
  of fighting it.
- Then item 4: the assembly inequality, one page, five certified constants, Lean socket.

*Drafted 2026-08-07 ~01:00 (session 2, after §802). Constants marked TBD are not certified;
nothing here banks until each B-piece has its arb line. Provenance: §791 envelope, §796/§800-§802
pencils, §798 KRD leads, §801 p⊕r form.*
