import numpy as np

# --- exact von Mangoldt over a window via linear sieve, delta = Lambda - 1 ---
def build_dLambda(N):
    Lam = np.zeros(N+1)
    spf = np.zeros(N+1, dtype=np.int64)  # smallest prime factor
    for i in range(2, N+1):
        if spf[i]==0:
            for j in range(i, N+1, i):
                if spf[j]==0: spf[j]=i
    for n in range(2, N+1):
        p = spf[n]; m = n
        # n is prime power iff repeatedly dividing by p reaches 1
        while m % p == 0: m //= p
        if m == 1:
            Lam[n] = np.log(p)
    return Lam - 1.0   # delta-Lambda, mean ~0

N = 2_000_000
dL = build_dLambda(N)
print(f"N={N}  E[dL]={dL[2:].mean():+.5f}  E[dL^2]={ (dL[2:]**2).mean():.4f}")

# accessible-odd triples: distinct positive offsets summing to s, all coprime-admissible mod 3,5,7
def admissible(t):
    for p in (3,5,7):
        if len(set(x % p for x in t)) == p:  # covers all residues -> non-admissible (a p always divides)
            return False
    return True
def triples_at_sum(s):
    out=[]
    for a in range(1,s):
        for b in range(a+1,s):
            c=s-a-b
            if c>b and admissible((a,b,c)): out.append((a,b,c))
    return out

def u_vec(t, base, span):
    n = np.arange(base, base+span)
    return dL[n+t[0]]*dL[n+t[1]]*dL[n+t[2]]

def gram(s, base=2, span=None):
    if span is None: span = N - 60 - base
    T = triples_at_sum(s)
    U = np.stack([u_vec(t,base,span) for t in T])   # rows = triples
    G = (U @ U.T)/span                               # sample Gram, PSD exactly
    d = np.sqrt(np.diag(G))
    R = G/np.outer(d,d)                              # unit-diagonal correlation
    ev = np.linalg.eigvalsh(R)
    off = R - np.eye(len(T))
    sigma_off = np.sqrt((off**2).sum()/(len(T)*(len(T)-1)))
    # least leading principal minor (Cholesky-equivalent positivity witness) of G
    try:
        L = np.linalg.cholesky(G); chol_ok=True; minpiv=float(np.min(np.diag(L))**2)
    except np.linalg.LinAlgError:
        chol_ok=False; minpiv=float('nan')
    return T, G, R, ev, sigma_off, chol_ok, minpiv, d

for s in (24,30,36,42):
    T,G,R,ev,so,ok,mp,d = gram(s)
    print(f"s={s:2d} n={len(T):3d}  lam_min(R)={ev.min():.4f}  lam_max={ev.max():.2f}  "
          f"sigma_off={so:.4f}  sig*sqrt(n)={so*np.sqrt(len(T)):.3f}  chol_PD={ok}  minpiv={mp:.3e}  d_max={d.max():.2f}")

print("\n=== N-stability of sigma_off (population vs finite-sample) ===")
for span_frac in (0.25, 0.5, 1.0):
    sp = int((N-60)*span_frac)
    T,G,R,ev,so,ok,mp,d = gram(s=36, base=2, span=sp)
    print(f"  span={sp:>9d}  lam_min(R)={ev.min():.4f}  sigma_off={so:.4f}")

print("\n=== matched Gaussian null: same (n,N), independent columns ===")
rng = np.random.default_rng(0)
for s in (30,36,42,48):
    T = triples_at_sum(s); n=len(T); span=N-60
    # ours
    U = np.stack([u_vec(t,2,span) for t in T]); G=(U@U.T)/span
    dd=np.sqrt(np.diag(G)); R=G/np.outer(dd,dd); lam_ours=np.linalg.eigvalsh(R).min()
    # matched null: n independent standard-normal streams, same sample count, sample-correlation edge
    Xn = rng.standard_normal((n, min(span,200000)))
    Rn = np.corrcoef(Xn); lam_null=np.linalg.eigvalsh(Rn).min()
    off=R-np.eye(n); so=np.sqrt((off**2).sum()/(n*(n-1)))
    print(f"  s={s:2d} n={n:3d}  lam_min_OURS={lam_ours:.4f}  lam_min_NULL(iid,200k)={lam_null:.4f}  sigma_off={so:.4f}  1-2*so*sqrt(n)={1-2*so*np.sqrt(n):.3f}")
