#!/usr/bin/env python3
"""§419: two-instrument cross-check of §417 (the discipline that caught §386). Replicate the cancellation-
domination finding on the INDEPENDENT Fejér/sinc² instrument (l2_fejer). If ||prime||_op/archmargin ALSO blows
up while push/||op||→0, the '292x energy hopeless, cancellation-dominated' crux is instrument-INDEPENDENT (real).
If Gaussian-specific, §417 is weakened."""
import numpy as np, l2_fejer as LF
BAND=LF.BAND
print(f"§419 Fejér-instrument cancellation check (cross-check of §417). band{BAND}.")
print(f"{'s':>5} {'n':>3} {'archmargin':>11} {'prime_push':>11} {'||prime||op':>11} {'push/||op||':>11} {'||op||/arch':>11}")
for s in np.arange(6.0,2.4,-0.3):
    cs=LF.probes(BAND,s)
    if len(cs)<4: continue
    f=LF.build(cs); AP=f['archpole']; PR=f['prime']
    la=float(np.linalg.eigvalsh(AP)[0])
    push=la-float(np.linalg.eigvalsh(AP-PR)[0])
    ev=np.linalg.eigvalsh(PR); opn=max(abs(ev[0]),abs(ev[-1]))
    print(f"{s:5.1f} {len(cs):3d} {la:+11.4e} {push:+11.4e} {opn:11.4e} {push/opn:11.3f} {opn/la if la>1e-15 else np.nan:11.3f}",flush=True)
print("\nreading: ||op||/arch blows up + push/||op||→0 on Fejér too => §417 cancellation-domination is instrument-")
print("independent (real RH-content, not a Gaussian artifact). Same pattern = the two-instrument rule confirms it.")
