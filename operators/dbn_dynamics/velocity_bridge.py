# Close the bridge BOTH-SIDED: log-gas velocity dx_j/ds = 2 sum_{k!=j} 1/(g_j-g_k)  ==  Xi''(g_j)/Xi'(g_j)
# (Hadamard). LHS from ZEROS (log-gas side), RHS from Xi=xi(1/2+iz) derivatives (PRIME side). Must match.
import mpmath as mp
mp.mp.dps=30
def xi(s): return mp.mpf('0.5')*s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
def Xi(z): return xi(mp.mpf('0.5')+1j*z)   # real for real z
# zeros
G=[]
with open('/Users/scones/projects/riemann/data/zeros_100k.txt') as f:
    for line in f:
        s=line.strip()
        if s:
            try: G.append(mp.mpf(s.split()[0]))
            except: pass
        if len(G)>=2000: break
G=G[:2000]
allz=[-g for g in G[::-1]]+list(G)   # +/- symmetric
print("j   gamma_j     LHS: 2*sum 1/(g_j-g_k) [zeros]     RHS: Xi''/Xi' [prime-side]      match")
for j in [0,1,2,3,4]:
    gj=G[j]
    # LHS: sum over all other zeros (+/-), truncated at 2000 each side (conditionally convergent, symmetric)
    lhs=mp.mpf(0)
    for gk in allz:
        if abs(gk-gj)>1e-20: lhs+=1/(gj-gk)
    lhs*=2
    # RHS: Xi''(gj)/Xi'(gj) via mpmath diff (real)
    Xp=mp.re(mp.diff(Xi,gj,1)); Xpp=mp.re(mp.diff(Xi,gj,2))
    rhs=Xpp/Xp
    print(f"{j+1}  {float(gj):.5f}   {mp.nstr(lhs,10):>22}   {mp.nstr(rhs,10):>22}   ratio={float(lhs/rhs):.4f}", flush=True)
print("\nmatch => log-gas velocity = Xi''/Xi' EXACT => dBN flow is PRIME-DETERMINED (bridge fully closed both-sided).")
