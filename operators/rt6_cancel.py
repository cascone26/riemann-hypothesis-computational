#!/usr/bin/env python3
"""§417A: how much does RH rely on CANCELLATION vs energy? The prime form can push λmin down by at most its
operator norm ||prime||_op (worst-case alignment). The ACTUAL push along the critical direction is prime_push.
ratio = prime_push/||prime||_op near 1 => primes align with the critical direction (no cancellation, tight);
ratio << 1 => the critical eigenvector is nearly ORTHOGONAL to the prime form's dominant modes (RH exploits
cancellation). Also ||prime||_op vs archmargin: if ||prime||_op > archmargin, the worst-case COULD violate,
so RH genuinely needs the cancellation (not a trivial size bound)."""
import numpy as np, l2_edge_sweep as L
BAND=(14.0,44.0)
print(f"§417A cancellation ratio. band{BAND} w=3.")
print(f"{'s':>5} {'n':>3} {'archmargin':>11} {'prime_push':>11} {'||prime||op':>11} {'push/||op||':>11} {'||op||/arch':>11}")
for s in np.arange(6.0,2.9,-0.3):
    cs=L.probes(BAND,s)
    if len(cs)<4: continue
    f=L.build_forms(cs,3.0); AP=f['archpole']; PR=f['prime']
    la=float(np.linalg.eigvalsh(AP)[0])
    push=la-float(np.linalg.eigvalsh(AP-PR)[0])
    opn=float(np.linalg.eigvalsh(PR)[-1])   # largest eigenvalue of prime (PSD? check via abs)
    opn=max(abs(np.linalg.eigvalsh(PR)[0]),abs(np.linalg.eigvalsh(PR)[-1]))
    print(f"{s:5.1f} {len(cs):3d} {la:+11.4e} {push:+11.4e} {opn:11.4e} {push/opn:11.3f} {opn/la:11.3f}",flush=True)
print("\nreading: push/||op|| → small => cancellation-dominated (critical dir avoids prime modes); →1 => aligned/tight.")
print("||op||/arch >1 => worst-case prime COULD exceed archmargin, so RH needs the cancellation genuinely.")
