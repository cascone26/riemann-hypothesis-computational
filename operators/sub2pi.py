#!/usr/bin/env python3
"""ANGLE 1 — IS THERE A PROVABLE POSITIVITY FOR 4W < 2π?
Sweep 4W ∈ {2, 3, 4, 5, 6, 6.5} across multiple windows/heights and characterize the corner.
Key questions:
  (1) Does a corner (#neg(M_σ)=1) exist for each 4W?
  (2) Does the negative eigenvalue depth μ change monotonically with 4W?
  (3) Is there a sharp transition at 4W=2π≈6.28 (does it vanish below?)?
  (4) If μ vanishes for small 4W, is M_σ⪰0 (no corner at all)?
Report actual #neg vs 4W data and whether 2π is a threshold."""

import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
import sys

def vonmangoldt(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam

NMAX=5000; LAM=vonmangoldt(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)

def analyze_corner(a, b, fourW, frac=0.65, dr=0.004):
    """Compute M_σ and characterize its spectrum, especially negative eigenvalues."""
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RGRID=np.arange(-(b+60),(b+60),dr)
    OMEGA=np.real(digamma(0.25+0.5j*RGRID))-np.log(np.pi); CEN=np.linspace(a,b,N)

    def fej(x):
        y=W*x; out=np.ones_like(y); nz=np.abs(y)>1e-12; out[nz]=(np.sin(y[nz])/y[nz])**2; return out

    P=np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G=(P@P.T)*dr
    arch=(P*OMEGA[None,:])@P.T*dr/(2*np.pi)

    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]; cosu=np.cos(np.outer(ln,RGRID)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cosu.T)*dr/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)

    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))

    # M_σ = compression of symbol σ(r) = Ω/2π - P(r)
    Pw=(1/np.pi)*(lw[:,None]*np.cos(np.outer(ln,RGRID))).sum(axis=0); sig=OMEGA/(2*np.pi)-Pw
    Msig=(P*sig[None,:])@P.T*dr; Msig=0.5*(Msig+Msig.T)

    # Whiten by G
    gw,gv=eigh(G); B=gv@np.diag(1/np.sqrt(gw))@gv.T
    Msw=B@Msig@B

    # Eigenvalues of M_σ
    ev_Msig=eigh(Msw,eigvals_only=True)
    nneg=int((ev_Msig<0).sum())

    if nneg > 0:
        # Deepest negative eigenvalue
        mu=-float(ev_Msig[0])  # -μ is the eigenvalue, μ>0 is its depth
    else:
        mu=0.0

    # Symbol properties
    sig_min=sig.min()
    sig_mean=sig.mean()
    neg_frac=float((sig<0).mean())

    return {
        'N': N,
        'nneg': nneg,
        'mu': mu,
        'ev_min': float(ev_Msig[0]),
        'ev_spectrum': ev_Msig[:min(5,N)],
        'sig_min': sig_min,
        'sig_mean': sig_mean,
        'neg_frac': neg_frac,
    }

print("="*130)
print("ANGLE 1: sweep 4W ∈ {2, 3, 4, 5, 6, 6.5} for corner existence and μ behavior across heights.")
print("=" * 130)
print(f"{'4W':>5} {'window':>16} {'N':>3} {'#neg(M_σ)':>10} {'μ=|λ_min|':>12} {'min(σ)':>11} {'σ mean':>10} {'σ negfrac':>10}")
print("-" * 130)

fourW_values = [2.0, 3.0, 4.0, 5.0, 6.0, 6.5]
windows_by_height = [
    (14, 44),    # lower height
    (30, 60),    # mid-lower
    (60, 90),    # mid
    (120, 150),  # mid-upper
    (300, 330),  # high
]

results_by_fourW = {fw: [] for fw in fourW_values}

for fourW in fourW_values:
    print(f"\n4W = {fourW:.1f}:")
    for a, b in windows_by_height:
        try:
            res = analyze_corner(a, b, fourW, frac=0.65)
            results_by_fourW[fourW].append((a, b, res))
            print(f"{fourW:>5.1f} {f'[{a:>3},{b:>3}]':>16} {res['N']:>3} {res['nneg']:>10} {res['mu']:>12.4e} {res['sig_min']:>+11.5f} {res['sig_mean']:>10.5f} {res['neg_frac']:>10.4f}")
        except Exception as e:
            print(f"{fourW:>5.1f} {f'[{a:>3},{b:>3}]':>16} ERROR: {str(e)[:60]}")

print("\n" + "="*130)
print("ANALYSIS: corner existence and μ trends across 4W values")
print("="*130)

for fourW in fourW_values:
    data = results_by_fourW[fourW]
    if not data:
        continue

    nneg_list = [res['nneg'] for _, _, res in data]
    mu_list = [res['mu'] for _, _, res in data]
    sig_min_list = [res['sig_min'] for _, _, res in data]

    corner_exists = any(nneg == 1 for nneg in nneg_list)
    all_psd = all(nneg == 0 for nneg in nneg_list)
    max_nneg = max(nneg_list)
    max_mu = max(mu_list)
    min_sig = min(sig_min_list)

    print(f"\n4W = {fourW:.2f} (vs 2π={2*np.pi:.3f}):")
    print(f"  corner (#neg=1) found: {corner_exists} | all M_σ⪰0: {all_psd} | max #neg: {max_nneg}")
    print(f"  max μ (deepest neg eig): {max_mu:.4e} | min σ: {min_sig:.5f}")
    print(f"  #neg across windows: {nneg_list}")
    print(f"  μ across windows: {[f'{m:.3e}' for m in mu_list]}")

print("\n" + "="*130)
print(f"KEY: 2π ≈ {2*np.pi:.3f}. Is there a threshold below 2π where M_σ becomes PSD (μ→0, #neg→0)?")
print("="*130)

# Analyze threshold behavior
print("\nTHRESHOLD ANALYSIS:")
for fw_idx in range(len(fourW_values)-1):
    fw1, fw2 = fourW_values[fw_idx], fourW_values[fw_idx+1]
    data1 = results_by_fourW[fw1]
    data2 = results_by_fourW[fw2]

    if not data1 or not data2:
        continue

    # Use the deepest corner window [60,90] if available
    res1 = next((res for a, b, res in data1 if a==60 and b==90), data1[-1])
    res2 = next((res for a, b, res in data2 if a==60 and b==90), data2[-1])

    mu1, nneg1 = res1['mu'], res1['nneg']
    mu2, nneg2 = res2['mu'], res2['nneg']

    print(f"\n4W: {fw1:.1f} → {fw2:.1f}")
    print(f"  #neg: {nneg1} → {nneg2}  |  μ: {mu1:.3e} → {mu2:.3e}  (ratio: {mu2/mu1 if mu1>0 else 'inf'})")
    if mu1 > 0 and mu2 == 0:
        print(f"  ⚡ TRANSITION: corner VANISHES between 4W={fw1} and 4W={fw2}")
    elif mu1 > 0 and mu2 > 0 and mu2 < mu1:
        print(f"  → μ DECREASING with 4W")
    elif mu1 > 0 and mu2 > 0 and mu2 > mu1:
        print(f"  → μ INCREASING with 4W (unexpected)")

print("\n" + "="*130)
