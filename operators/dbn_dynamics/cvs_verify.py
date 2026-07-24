import sys, connes_cvs as m, mpmath as mp
c=13; N=int(sys.argv[1]); dps=int(sys.argv[2]); T=int(sys.argv[3]) if len(sys.argv)>3 else 400
mp.mp.dps=dps
Q=m.build_galerkin_matrix(c, N=N, T=T, dps=dps)
lam, vec = m.compute_ground_state(Q)
print(f"c={c} N={N} dps={dps} T={T}: lambda_min = {mp.nstr(lam,6)}  (Groskin target N=100: 2.865e-59)")
L=float(mp.log(c))
zs=m.extract_zeros(vec, L, n_zeros=6, dps=dps)
print("ground-state zeros:", [mp.nstr(z,8) for z in zs])
print("gamma_1..4 = 14.134725, 21.022040, 25.010858, 30.424876")
