# Verify cascade: eigenvector v_hat is O(1) OFF zeros, ~0 AT zeros (real signal, not uniformly tiny)
import mpmath as mp
mp.mp.dps=35
import sys; sys.path.insert(0,'/tmp')
exec(open('/tmp/cascade_test.py').read().split('print("CASCADE')[0])  # reuse builder up to eigendecomp
order=sorted(range(D),key=lambda k: Ev[k])
gam=[14.134725,21.022040,25.010858]; nong=[7.0,11.0,17.5,27.0,35.0]
for rank in [0,3,6]:
    k=order[rank]; c=[Q[i,k] for i in range(D)]
    atg=[float(vhat(c,mp.mpf(g))) for g in gam]
    atn=[float(vhat(c,mp.mpf(z))) for z in nong]
    print(f"rank {rank} (eig {mp.nstr(Ev[k],3)}): |vhat| AT zeros={[f'{x:.2e}' for x in atg]}  OFF zeros={[f'{x:.2e}' for x in atn]}",flush=True)
print("(AT<<OFF => genuine vanishing at zeros)",flush=True)
