import numpy as np
exec(open('weyl_wkb_tail_test.py').read().split('gs = np.linspace')[0])
Z = np.load('zeros500.npy')
zin = Z[(Z>=12)&(Z<=32)]
gs = np.linspace(12,32,60)
import sys
X = float(sys.argv[1])
for name, fn in [('naive', B_naive), ('wkb', B_wkb)]:
    r = roots_from(fn, X, gs)
    nn = min(len(r), len(zin), 4)
    if nn == 0:
        print(f"{X:5.0f} | {name:8} | {len(r):6d} | (no roots)"); continue
    errs = 100*(r[:nn]-zin[:nn])/zin[:nn]
    print(f"{X:5.0f} | {name:8} | {len(r):6d} | {np.mean(np.abs(errs)):8.3f}%   " +
          ", ".join(f"{e:+.2f}%" for e in errs))
