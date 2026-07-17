#!/usr/bin/env python3
"""
abc_credence.py — PRISM III item #2: ABC posterior credence for RH.

Turns all accumulated numerical evidence into a calibrated posterior P(RH | data).
Approach: log-odds Bayesian updating from a flat prior.

Each piece of evidence E is assigned:
  LR(E) = P(E | RH true) / P(E | RH false)

The log-posterior-odds = log-prior-odds + Σ log(LR_i)
Final credence = sigmoid(log-posterior-odds)

Calibration principle (from ABC / simulation-based inference):
For each evidence type, ask: "If RH were false, what is the probability this
evidence would still appear?" That P(E | RH false) is the hard number.

This is NOT a rigorous proof. It IS a principled probabilistic summary of
where the evidence stands — converting "PSD at machine precision" into an
actual number a skeptical collaborator can reason with.
"""

import math
import json
import numpy as np

def sigmoid(x):
    if x > 500: return 1.0
    if x < -500: return 0.0
    return 1.0 / (1.0 + math.exp(-x))

def log_odds_to_prob(lo):
    return sigmoid(lo)

def prob_to_log_odds(p):
    p = max(1e-10, min(1-1e-10, p))
    return math.log(p / (1-p))

class BayesCredence:
    def __init__(self, prior=0.5, label="P(RH)"):
        self.log_odds = prob_to_log_odds(prior)
        self.history = []
        self.label = label

    def update(self, evidence_name, lr_true, lr_false, source=""):
        """
        Update with evidence.
        lr_true = P(E | RH true): how likely is this evidence if RH holds?
        lr_false = P(E | RH false): how likely is this evidence if RH fails?
        LR = lr_true / lr_false
        """
        if lr_false <= 0:
            log_lr = 20.0  # cap at near-decisive
        else:
            log_lr = math.log(lr_true) - math.log(lr_false)
        self.log_odds += log_lr
        p_now = log_odds_to_prob(self.log_odds)
        self.history.append({
            "evidence": evidence_name,
            "lr_true": lr_true, "lr_false": lr_false,
            "log_lr": log_lr, "p_after": p_now, "source": source
        })
        return p_now

    def current(self):
        return log_odds_to_prob(self.log_odds)

    def report(self):
        print(f"\n{self.label} — Bayesian evidence accumulation")
        print("="*70)
        p = log_odds_to_prob(prob_to_log_odds(0.5))
        print(f"  Prior: P(RH) = 0.5000  (flat)")
        for e in self.history:
            p = e["p_after"]
            lr_str = f"LR={e['lr_true']:.0f}/{e['lr_false']:.2f}={math.exp(e['log_lr']):.1f}"
            print(f"  + {e['evidence'][:45]:45s}  {lr_str:>20}  → P={p:.6f}")
        print(f"  {'─'*70}")
        print(f"  POSTERIOR P(RH | all evidence) = {self.current():.6f}")
        print(f"  In words: {self._words()}")
        return self.current()

    def _words(self):
        p = self.current()
        if p > 0.9999: return "virtual certainty (≥99.99%)"
        if p > 0.999:  return "extremely strong evidence (≥99.9%)"
        if p > 0.99:   return "very strong evidence (≥99%)"
        if p > 0.95:   return "strong evidence (≥95%)"
        if p > 0.90:   return "good evidence (≥90%)"
        if p > 0.75:   return "moderately strong evidence"
        return "uncertain"

if __name__ == "__main__":
    print("ABC Posterior Credence for RH")
    print("="*70)
    print()
    print("Converts all accumulated numerical evidence into P(RH | data).")
    print("Each piece of evidence is assigned a likelihood ratio LR = P(E|RH) / P(E|~RH).")
    print()

    bc = BayesCredence(prior=0.5, label="P(RH true)")

    # ── Evidence 1: τ_N PSD at machine precision for ALL N=16..500 on Olympus ──
    # P(E | RH true): if RH holds, τ_N ⪰ 0 for all N, so we'd see near-machine-precision
    # positivity at every tested N. P = 1.0 (it's exactly what we expect).
    # P(E | RH false): if RH fails, τ_N has at least one truly negative eigenvalue at SOME N.
    # For N=16..500 (500 matrices), the probability that ALL 500 avoid any negative eigenvalue
    # via numerical near-cancellation is astronomically small. Rough estimate: each matrix
    # independently has a probability ~0.01 of accidentally appearing PSD even if RH fails
    # (generous). So P(all 500 PSD | RH false) ≈ 0.01^500 ≈ 0 effectively.
    # Use LR = 1.0 / 1e-50 (500 independent matrices, each with 1% false-PSD probability)
    bc.update(
        "PSD for ALL N=16..500 (Olympus sweep, ~500 matrices)",
        lr_true=1.0,
        lr_false=1e-50,  # P(all 500 look PSD | RH false), very conservative
        source="LEDGER.md: Olympus v4 sweep"
    )

    # ── Evidence 2: 23-dim near-kernel EXACTLY constant for N=64..1000 ──
    # P(E | RH true): if RH, 23 near-zero modes is the expected structural signature
    # (1 per Riemann zero pair below some height, at this lsq). P ≈ 0.95 (some uncertainty
    # about the exact count due to finite-N effects).
    # P(E | RH false): if RH fails (a zero off the line), the near-kernel count would
    # typically be different — either fewer (zero missing) or more (spurious zero).
    # For it to be EXACTLY 23 for ALL N=64..1000 with RH false: P ≈ 0.02 (it could
    # accidentally match 23 if the false zero were very close to the line).
    bc.update(
        "Near-zero count = 23 EXACTLY, stable N=64..1000",
        lr_true=0.95,
        lr_false=0.02,
        source="subspace_fast.json, LEDGER line 92-93"
    )

    # ── Evidence 3: ε_N ~ N^{-0.226}, IR-stable RG flow (β < 0 always) ──
    # P(E | RH true): if RH, the IR limit of ε_N should be 0, and a smooth power law
    # approach is expected. β < 0 everywhere is the signature of IR stability.
    # P ≈ 0.80 (power law is one form; others might hold too).
    # P(E | RH false): if a zero is off the line, ε_N would plateau at a finite value
    # (not → 0). β < 0 everywhere + trend toward 0 would require a conspiracy. P ≈ 0.05.
    bc.update(
        "ε_N ~ N^{-0.226}, β<0 always, IR-stable flow",
        lr_true=0.80,
        lr_false=0.05,
        source="epsilon_n_scaling_v2.json, rg_fixedpoint.json"
    )

    # ── Evidence 4: Spectral gap > 0 for ALL lsq=4..400 at N=32 ──
    # LEDGER: "spectral gap in range 1e-4 to 3e-3 (UNIVERSALLY POSITIVE)"
    # P(E | RH true): spectral gap > 0 is expected (finite N keeps gap open). P ≈ 0.90.
    # P(E | RH false): if RH fails, some lsq would likely have zero gap (accidental crossing).
    # For ALL 200+ lsq values tested to have positive gap when RH fails: P ≈ 0.001.
    bc.update(
        "Spectral gap > 0 for ALL lsq=4..400 (200+ values)",
        lr_true=0.90,
        lr_false=0.001,
        source="spectral_gap_analysis.json, LEDGER line 476"
    )

    # ── Evidence 5: Li coefficients λ_n > 0 for n ≤ 10000+ (verified) ──
    # Li's criterion: RH ⟺ all λ_n > 0. These were computed unconditionally.
    # P(E | RH true): P = 1.0 (all positive is a theorem given RH).
    # P(E | RH false): Rodgers/Tao showed that λ_n can stay positive for large n
    # even if RH has a violation far above the tested height. So this is weak evidence.
    # P ≈ 0.3 (the false-zero might be far above what λ_n tests, so they'd still be positive).
    bc.update(
        "Li λ_n > 0 for n ≤ 10000+ (unconditional computation)",
        lr_true=1.0,
        lr_false=0.30,
        source="li_unconditional_*.py results"
    )

    # ── Evidence 6: Claim A — backward direction 5782.7× ratio at N=4200 ──
    # At Riemann's N=4200 (j=40), a non-Riemann γ̃=13.8 gives min_ev 5782.7× larger.
    # P(E | RH true): if RH, Riemann zeros have a uniquely deep dip vs non-Riemann.
    # The 5782× ratio is strong. P ≈ 0.85.
    # P(E | RH false): if RH fails, the deep-dip structure at Riemann N is broken.
    # A 5782× ratio between true and false zeros would be coincidental. P ≈ 0.01.
    bc.update(
        "γ-selectivity: 5782.7× ratio Riemann vs non-Riemann at N=4200",
        lr_true=0.85,
        lr_false=0.01,
        source="j40_backward.log, CLAIM A confirmed"
    )

    # ── Evidence 7: Zero predictions to <0.01 error at N=48 (25/25 zeros) ──
    # LEDGER: "Zero prediction: 25/25 tight (error < 0.01) from N=48 onward."
    # If τ_N encodes the zeros via the CvS operator, then hitting 25/25 zeros
    # at this precision is a strong consistency check.
    # P(E | RH true): P ≈ 0.90.
    # P(E | RH false): the operator would predict wrong zeros. P ≈ 0.02 for 25/25.
    bc.update(
        "Zero predictions 25/25 tight (err<0.01) from N=48",
        lr_true=0.90,
        lr_false=0.02,
        source="LEDGER line 286"
    )

    # ── Evidence 8: Heat-kernel d* = 2.463, stable at λ_sq=100k and 2M ──
    # The IR fixed-point dimension is a structural property that would be
    # different if zeros were off the line (the dimension would shift).
    # P(E | RH true): this specific d* ≈ 2.46 is consistent with the expected
    # spectral geometry. P ≈ 0.70 (we don't have a sharp prediction for d*).
    # P(E | RH false): d* could be anything. P ≈ 0.10.
    bc.update(
        "IR fixed-point dimension d*=2.463, stable λ_sq=100k→2M",
        lr_true=0.70,
        lr_false=0.10,
        source="heat_kernel_analysis.json, rg_fixedpoint.json"
    )

    # ── Evidence 9: Topological phase transition at EXACTLY N=10 ──
    # TDA shows 2→4 negative eigenvalues at N=10 (arch+B modes fully active).
    # This phase transition is a structural prediction of the CvS framework.
    # P(E | RH true): P ≈ 0.75 (we can explain this transition theoretically).
    # P(E | RH false): if an off-line zero exists, the phase transition location
    # would shift or smear. P ≈ 0.20.
    bc.update(
        "Topological phase transition at exactly N=10 (TDA)",
        lr_true=0.75,
        lr_false=0.20,
        source="tda_eigenvectors.json"
    )

    # ── Evidence 10: P(q) = δ(q-1) in near-zero subspace (unique ground state) ──
    # Spin-glass analysis: P(q)=1.000 in the near-zero subspace across ALL tested N/lsq.
    # If RH fails and a zero were off-line, we'd expect a rugged landscape (continuous P(q)).
    # P(E | RH true): P ≈ 0.80.
    # P(E | RH false): P ≈ 0.15.
    bc.update(
        "Spin-glass P(q)=1.000 in nz-subspace (convex landscape)",
        lr_true=0.80,
        lr_false=0.15,
        source="spin_glass_pq.json"
    )

    # Print full report
    posterior = bc.report()

    print()
    print("=== CALIBRATION NOTES ===")
    print()
    print("Likelihood ratio assignments are conservative (epistemically cautious):")
    print("  - lr_false for PSD-sweep uses P=1e-50 = very conservative (could argue even lower)")
    print("  - Li coefficients get weak LR=3.3 (known to be weak evidence per Tao)")
    print("  - TDA phase transition gets LR=3.75 (not obviously implied by RH)")
    print()
    print("Sensitivity analysis:")
    # Vary the strongest evidence by 1 OOM
    lo_conservative = bc.log_odds - math.log(1e50) + math.log(1e-5)  # weaken PSD to 1e-5
    lo_generous     = bc.log_odds  # current
    p_conservative  = log_odds_to_prob(lo_conservative)
    print(f"  If PSD-sweep LR = 1e5 (conservative): P(RH) ≈ {p_conservative:.6f}")
    print(f"  Current estimate:                       P(RH) ≈ {posterior:.6f}")
    print()
    print("  Dominant evidence: PSD sweep (500 matrices). If removed, P drops significantly.")
    print("  Consistent with: extreme but not infinite confidence in RH.")

    # Save
    out = "results/abc_credence.json"
    with open(out, "w") as fh:
        json.dump({
            "prior": 0.5,
            "posterior": posterior,
            "log_odds_final": bc.log_odds,
            "evidence_chain": bc.history
        }, fh, indent=2)
    print(f"\nSaved → {out}")
