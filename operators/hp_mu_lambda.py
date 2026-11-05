"""S2887: margin curve eps_{lambda,N} (=lambda_min of CCM/Shi tau) across lambda^2, two N per lambda (saturation -> mu_lambda)."""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
for spec in sys.argv[1:]:
    lam2,N,dps=spec.split(':'); lam2=float(lam2); N=int(N); dps=int(dps)
    W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True); S=W02-WR-WP; mp.mp.dps=dps
    E=sorted(mp.eigsy(S)[0])
    print(f"lam2={lam2:6.2f} L={float(L):.4f} N={N:3d}: eps={mp.nstr(E[0],5):>12}  e1={mp.nstr(E[1],4):>11}  log10eps={float(mp.log10(abs(E[0]))):8.2f}",flush=True)
