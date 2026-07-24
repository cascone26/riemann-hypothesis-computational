#!/usr/bin/env python3
"""§176: ANALYTICAL prime-gate — per-prime Hellmann-Feynman decomposition of the near-null.
From Eq 2.11 [25]: Q = arch − Σ_{n≤c}Λ(n)⟨f|V(n)f⟩. First-order gate δλ_min(n) = ⟨ψ|V_n|ψ⟩ where V_n = the
matrix built from prime-power n's contribution alone, ψ = near-null eigenvector. Eigenvector expectation ⟹
ROBUST (§173: eigenvector stable to 1e-100, unlike the tiny eigenvalue — sidesteps the §170 noise wall).
Reports h(n)=⟨ψ|V_n|ψ⟩ per prime power, vs von Mangoldt weight Λ(n)/√n, sign, and Σ (reconstruct total lift)."""
import sys, mpmath as mp
sys.path.insert(0, ".venv_cvs/lib/python3.14/site-packages")
import connes_cvs.operator as O
from cvs_fast import precompute_components

c   = int(sys.argv[1]) if len(sys.argv)>1 else 13
N   = int(sys.argv[2]) if len(sys.argv)>2 else 80
dps = int(sys.argv[3]) if len(sys.argv)>3 else 150
T   = int(sys.argv[4]) if len(sys.argv)>4 else 400
mp.mp.dps=dps
print(f"§176 prime-gate H-F decomposition — c={c} N={N} dps={dps} T={T}\n", flush=True)

comp = precompute_components(c, N, T, dps)
idx = comp["idx"]; DIM = 2*N+1
# even-sector projector
V = mp.matrix(DIM, N+1); V[N,0]=mp.mpf(1); r2=1/mp.sqrt(2)
for k in range(1,N+1): V[N+k,k]=r2; V[N-k,k]=r2

def assemble_from(valmap, dermap):
    """even-sector matrix from a psi value/deriv map (no arch/pole if only prime given)."""
    Q=mp.matrix(DIM,DIM)
    for i in range(DIM):
        m=i-N
        for j in range(DIM):
            nn=j-N
            Q[i,j]= dermap[nn] if m==nn else (valmap[m]-valmap[nn])/(m-nn)
    for i in range(DIM):
        for j in range(i+1,DIM):
            a=(Q[i,j]+Q[j,i])/2; Q[i,j]=a; Q[j,i]=a
    return V.T*Q*V

# full operator (arch+pole+all primes) -> near-null eigenvector psi
keep=[nn for nn in comp["pp"]]
psi_full={n: sum((comp["pp"][nn][0][n] for nn in keep), mp.mpf(0)) + comp["arcpol"][n] for n in idx}
psid_full={n: sum((comp["pp"][nn][1][n] for nn in keep), mp.mpf(0)) + comp["arcpol_d"][n] for n in idx}
Qe = assemble_from(psi_full, psid_full)
eigs, vecs = mp.eigsy(Qe)
k0 = min(range(N+1), key=lambda i: eigs[i])
psi = mp.matrix([vecs[r,k0] for r in range(N+1)])
lam0 = eigs[k0]
print(f"near-null λ_min = {mp.nstr(lam0,6)} (log10={float(mp.log10(abs(lam0))):.2f}); ψ = ground even eigenvector\n")

# per-prime-power Hellmann-Feynman contribution h(n) = <ψ| V_n |ψ>
def vmangoldt(n):
    # Λ(n)/√n : n=p^k -> log p / sqrt(n)
    p=2
    while n%p!=0: p+=1
    return float(mp.log(p)/mp.sqrt(n))
print(f"{'n(pp)':>6} {'Λ(n)/√n':>10} {'h(n)=<ψ|V_n|ψ>':>18} {'sign':>5}")
tot=mp.mpf(0)
rows=[]
for nn in sorted(comp["pp"]):
    Vn = assemble_from({m: comp["pp"][nn][0][m] for m in idx}, {m: comp["pp"][nn][1][m] for m in idx})
    hn = (psi.T*(Vn*psi))[0]
    tot+=hn
    rows.append((nn,hn))
    print(f"{nn:>6} {vmangoldt(nn):>10.4f} {mp.nstr(hn,6):>18} {'+' if hn>0 else '-':>5}", flush=True)
print(f"\nΣ h(n) over prime powers = {mp.nstr(tot,6)}  (total first-order prime contribution to λ_min on ψ)")
arch_exp = (psi.T*(assemble_from(comp['arcpol'],comp['arcpol_d'])*psi))[0]
print(f"<ψ|arch+pole|ψ> = {mp.nstr(arch_exp,6)};  arch + Σh = {mp.nstr(arch_exp+tot,6)}  (should ≈ λ_min = {mp.nstr(lam0,6)})")
print("\nFULL-PRECISION TELESCOPING CASCADE — residual R_m = arch + Σ_{n≤m} h(n) after each prime power:")
R = arch_exp
print(f"   {'arch':>6}: log10|R| = {float(mp.log10(abs(R))):>8.3f}")
for nn, hn in rows:
    R = R + hn
    if R != 0:
        print(f"   +n={nn:>2}: log10|R| = {float(mp.log10(abs(R))):>8.3f}   (|h(n)|=10^{float(mp.log10(abs(hn))) if hn!=0 else 0:.2f})")
print("If log10|R| drops monotonically ~linearly toward log10(λ_min)≈−58.5, each prime power peels off another factor:")
print("the near-null e^{−4πc} is a CASCADE of ~π(c) prime cancellations. That is the analytical anatomy.")
print("\nREAD: h(n) is the per-prime GATE (first-order δλ from prime n). Structure vs Λ(n)/√n, signs, and how Σ builds")
print("the near-null from the negative arch expectation = the analytical anatomy of Weil positivity. Robust (eigenvector).")
