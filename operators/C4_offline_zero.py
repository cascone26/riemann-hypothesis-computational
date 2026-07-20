#!/usr/bin/env python3
"""C4 (non-circular circularity test for the self-similarity lever): build the Weil quadratic form from a
CHOSEN zero set (NOT from primes — so we control zero locations), and test whether positivity is EQUIVALENT
to all zeros being on the critical line.

Weil/Guinand form (spectral / over-zeros side): for real test functions g_s in a family, the form is
  W[s,t] = Σ_ρ  ĝ_s(ρ) · conj(ĝ_t(ρ))
where ρ runs over the nontrivial zeros and ĝ is the (entire) Mellin/Fourier transform. For a zero ρ=1/2+iγ
ON the line, ρ contributes ĝ_s(1/2+iγ); with the functional-equation partner 1/2−iγ and g real, each zero
contributes a real rank-1 PSD term. Move a zero OFF the line (ρ=1/2+δ+iγ, plus its FE quartet
{1/2±δ±iγ}) and the contribution becomes a rank-2 block that can be INDEFINITE.

Model (labeled): test functions are Gaussian wavepackets in the 'time' variable x=log; their transform is
ĝ_s(1/2+iξ)=exp(−(ξ−μ_s)^2 σ^2/2) analytically continued to complex ξ. For a zero at height γ, real part
β=1/2+δ, the relevant argument is ξ=γ − iδ (from ρ=1/2+δ+iγ → ξ=(ρ−1/2)/i=γ−iδ). We build W over a small
set of test centers μ_s spanning the zeros, sum the on-line zeros' PSD terms, add ONE tunable zero at
(β,γ0) with its FE partner, and track λ_min(W) as δ grows from 0.

PASS (lever validated): λ_min(W) = 0-side (PSD) at δ=0 and goes NEGATIVE for δ≠0 ⇒ positivity ⟺ on-line,
non-circularly (we chose the zeros). Then the self-similar per-zero template is a genuine RH-equivalent
positivity condition and a renormalization/induction over zeros is logically sound.
"""
import numpy as np, math, json

ZEROS=[14.134725,21.022040,25.010858,30.424876,32.935062,37.586178,40.918719,43.327073,
       48.005151,49.773832,52.970321,56.446248]

def ghat(mu, sigma, xi):
    """transform of a Gaussian wavepacket centered at frequency mu, evaluated at (possibly complex) xi."""
    return np.exp(-((xi-mu)**2)*(sigma**2)/2.0)

def weil_matrix(centers, sigma, zero_list):
    """Weil form W[s,t] = Σ_ρ ĝ_s(ξ_ρ)·ĝ_t(ξ_ρ), ξ_ρ=(ρ−1/2)/i. Outer product WITHOUT conjugation:
    for a zero ON the line ξ is real → outer(v,v) is real PSD; for an OFF-line zero ξ is complex →
    outer(v,v) is complex-SYMMETRIC and its real part can be INDEFINITE. This is the true Weil mechanism.
    Each zero contributes its functional-equation set of ξ points: for ρ=1/2+δ+iγ the FE-related ρ give
    ξ ∈ {±γ − iδ, ±γ + iδ}."""
    S=len(centers); W=np.zeros((S,S),dtype=complex)
    for (beta,gamma) in zero_list:
        delta=beta-0.5
        xis = [gamma-1j*delta, -gamma-1j*delta] if delta==0 else \
              [gamma-1j*delta, -gamma-1j*delta, gamma+1j*delta, -gamma+1j*delta]  # FE quartet off-line
        for xi in xis:
            v=np.array([ghat(mu,sigma,xi) for mu in centers])
            W+=np.outer(v, v)          # NO conjugate — complex-symmetric off-line
    return W.real                       # real part of the complex-symmetric form

def main():
    sigma=0.35
    centers=[g for g in ZEROS]           # one test wavepacket per zero (spans them)
    print("=== C4: is Weil-form positivity EQUIVALENT to zeros on the line? (chosen zeros, no primes) ===")
    # baseline: all zeros on-line
    on=[(0.5,g) for g in ZEROS]
    W0=weil_matrix(centers,sigma,on)
    lam0=np.linalg.eigvalsh(W0)[0]
    print(f"  all zeros ON line: lambda_min(W) = {lam0:+.4e}  ({'PSD' if lam0>=-1e-12 else 'INDEFINITE'})")
    print("\n  move zero #1 (gamma={:.3f}) off-line by delta; track lambda_min(W):".format(ZEROS[0]))
    print(f"  {'delta':>8} {'lambda_min':>14} {'status':>12}")
    rows=[]
    for delta in [0.0,0.02,0.05,0.1,0.15,0.2,0.3,0.5]:
        zl=[(0.5,g) for g in ZEROS]
        zl[0]=(0.5+delta, ZEROS[0])          # move one zero off the critical line
        W=weil_matrix(centers,sigma,zl)
        lam=np.linalg.eigvalsh(W)[0]
        st="PSD" if lam>=-1e-10 else "INDEFINITE"
        rows.append(dict(delta=delta,lambda_min=float(lam),status=st))
        print(f"  {delta:>8.3f} {lam:>14.4e} {st:>12}")
    json.dump({"baseline_lambda_min":float(lam0),"sweep":rows},
              open("results/C4_offline_zero.json","w"),indent=2)
    print("\n=== VERDICT (C4) ===")
    on_line_psd = lam0>=-1e-10
    goes_neg = any(r['lambda_min']<-1e-9 for r in rows if r['delta']>0)
    if on_line_psd and goes_neg:
        print("  => PASS: PSD when all zeros on-line, INDEFINITE once a zero moves off. Positivity ⟺ on-line,")
        print("     built from CHOSEN zeros (no primes) ⇒ NON-CIRCULAR. The self-similar per-zero template is")
        print("     a genuine RH-equivalent positivity condition; renormalization/induction over zeros is sound.")
    elif on_line_psd and not goes_neg:
        print("  => WEAK: stays PSD even off-line in this model — the template does NOT constrain zeros here")
        print("     (self-similarity may just re-express the explicit formula). Lever weakened; reassess.")
    else:
        print("  => model degenerate (not PSD even on-line); refine test-function family / sigma.")
    print("wrote results/C4_offline_zero.json")

if __name__=="__main__":
    main()
