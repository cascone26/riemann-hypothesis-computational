#!/usr/bin/env python3
"""§429b VERIFY: is the ρ_edge→1 trend real, or a ZERO-TRUNCATION artifact (instrument has only 29 zeros, γ<101,
so high bands miss the zeros above them)? Check the explicit-formula residual ||full − zero||_max per band — the
correctness certificate. If residual is small (~1e-8) for all bands, ρ_edge trend is real. If it BLOWS UP for high
bands, the instrument is invalid there and ρ_edge<... is truncation garbage."""
import numpy as np, l2_fejer as LF
ZER=LF.ZEROS
def nz(band): return int(np.sum((ZER>=band[0])&(ZER<=band[1])))
print(f"§429b explicit-formula residual per band (correctness gate). zeros up to {ZER.max():.1f}, count={len(ZER)}")
print(f"{'band':>14} {'Nz_in':>5} {'Nz_above':>8} {'n':>3} {'||full−zero||':>13} {'archmargin':>11} {'ρ_edge':>7} {'VALID?':>7}")
for band in [(14,44),(20,50),(30,60),(40,70),(50,80),(60,90),(70,99)]:
    N=nz(band)
    if N<3: continue
    nz_above=int(np.sum(ZER>band[1]))
    for s in np.linspace((band[1]-band[0])/(N+3),(band[1]-band[0])/max(1,(N-1)),40):
        cs=LF.probes(band,s)
        if len(cs)==N:
            f=LF.build(cs)
            resid=np.max(np.abs(f['full']-f['zero']))
            am=float(np.linalg.eigvalsh(f['archpole'])[0])
            ev=np.linalg.eigvalsh(f['prime']); opn=max(abs(ev[0]),abs(ev[-1]))
            valid = resid<1e-6
            print(f"{str(band):>14} {N:5d} {nz_above:8d} {len(cs):3d} {resid:13.2e} {am:+11.4e} {opn/am:7.3f} "
                  f"{'YES' if valid else 'NO':>7}",flush=True)
            break
print("""
reading: residual small (VALID) for all bands => ρ_edge trend is real quantitative content. residual BLOWS UP
for high bands (few zeros above) => the explicit formula doesn't close there (truncation); the high-band ρ_edge<1
is an ARTIFACT of missing zeros, NOT a real unconditional zone. The valid range is where Nz_above is large enough.
""")
