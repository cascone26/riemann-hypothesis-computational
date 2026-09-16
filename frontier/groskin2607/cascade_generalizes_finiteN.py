"""
riemann11 -- test whether the §176/176b telescoping prime-cancellation cascade
(discovered on the LARGE-N near-null lambda_a operator) also governs the
FINITE-(c,N) Groskin/CvS margin -- an attempt at an ANALYTIC mechanism for the
empirical margin(c,N) ~ K*c^(-0.44N) law (margin_law.py/margin_law_results.txt),
using the exact same tau = W02 - WR - Wp assembly as arb_ldlt_certify.py, but at
float64 precision (REPORT-tier exploration, not Arb-rigorous).
"""
import numpy as np
import mpmath as mp

def prime_powers_up_to(c):
    primes = []
    x = 2
    while x <= c:
        if all(x % p for p in primes):
            primes.append(x)
        x += 1
    out = []
    for p in primes:
        q = p
        while q <= c:
            out.append((q, p))
            q *= p
    return out

def build_tau_and_terms(c, N):
    """Return (W02-WR block, list of per-prime-power Wp_idx matrices, pdata)."""
    L = np.log(c)
    PI = np.pi
    dim = 2*N+1
    ns = np.arange(-N, N+1)

    # S, CC, XC via mpmath digamma/trigamma (float-precision, small N so fine)
    S = np.zeros(N+1); CC = np.zeros(N+1); XC = np.zeros(N+1)
    for n in range(N+1):
        w = 2*PI*n/L
        z = complex(0.25, PI*n/L)
        psi = complex(mp.digamma(z))
        psi1 = complex(mp.polygamma(1, z))
        # geometric sums (finite truncation, k up to say 200 terms -> converges fast, e^{-c_k L})
        gS=gCC=gX1=gX2=0.0
        k=0
        while True:
            ck = 2*k+0.5
            e = np.exp(-ck*L)
            den = ck*ck + w*w
            gS += e/den
            if n!=0: gCC += e*w*w/(ck*den)
            gX1 += e*ck/den
            gX2 += e*(ck*ck-w*w)/(den*den)
            if e < 1e-18 and k>2: break
            k+=1
            if k>500: break
        if n==0:
            S[n]=0.0; CC[n]=0.0
        else:
            S[n] = 0.5*psi.imag - w*gS
            CC[n] = -0.5*(psi.real - float(mp.digamma(0.25))) + gCC
        XC[n] = 0.25*psi1.real - L*gX1 - gX2

    def Ssigned(nn):
        return S[nn] if nn>=0 else -S[-nn]

    sp2 = 16*PI*PI
    l2 = L*L
    pref02 = 32*L*np.sinh(L/4)**2
    U = np.exp(L/2)
    J = -2*np.log(U+1) + np.log(U*U+1) + 2*np.arctan(U) + np.log(2) - PI/2
    eL = np.exp(L)
    kappa = np.log(4*PI*(eL-1)/(eL+1)) + float(mp.euler)

    W02_WR = np.zeros((dim,dim))
    for i in range(dim):
        n = i-N
        for j in range(i,dim):
            m = j-N
            num = l2 - sp2*m*n
            den = (l2+sp2*m*m)*(l2+sp2*n*n)
            W02 = pref02*num/den
            if n==m:
                WR = kappa + 2*CC[abs(n)] + J - (2/L)*XC[abs(n)]
            else:
                WR = (Ssigned(m)-Ssigned(n))/(PI*(n-m))
            val = W02-WR
            W02_WR[i,j]=val; W02_WR[j,i]=val

    pdata = prime_powers_up_to(c)
    Wp_terms = []
    for (q,p) in pdata:
        y = np.log(q)
        weight = np.log(p)/np.sqrt(q)
        mat = np.zeros((dim,dim))
        for i in range(dim):
            n = i-N
            for j in range(i,dim):
                m = j-N
                if n==m:
                    val = 2*(1-y/L)*np.cos(2*PI*n*y/L)
                else:
                    val = (np.sin(2*PI*m*y/L)-np.sin(2*PI*n*y/L))/(PI*(n-m))
                mat[i,j]=val; mat[j,i]=val
        Wp_terms.append(weight*mat)
    return W02_WR, Wp_terms, pdata

def run(c, N):
    W02_WR, Wp_terms, pdata = build_tau_and_terms(c, N)
    Wp_total = sum(Wp_terms)
    tau = W02_WR - Wp_total
    evals, evecs = np.linalg.eigh(tau)
    idx = np.argmin(np.abs(evals))
    lam = evals[idx]
    psi = evecs[:, idx]
    psi = psi/np.linalg.norm(psi)
    h_arch = psi @ W02_WR @ psi
    print(f"c={c} N={N} dim={2*N+1}  lambda_min={lam:.6e} (log10={np.log10(abs(lam)):.2f})")
    print(f"  <psi|W02-WR|psi> = {h_arch:.6f}")
    running = h_arch
    print(f"  {'prime^k':>8} {'y=logq':>8} {'h(q)':>14} {'residual_after':>16} {'log10|resid|':>13}")
    for (q,p), mat in zip(pdata, Wp_terms):
        h = -(psi @ mat @ psi)   # tau has -Wp, so each prime's contribution to tau's expectation is -Wp_idx
        running += h
        print(f"  {q:>8} {np.log(q):>8.3f} {h:>14.6e} {running:>16.6e} {np.log10(abs(running)+1e-300):>13.2f}")
    print(f"  FINAL running residual = {running:.6e}  vs true lambda_min = {lam:.6e}  (ratio {running/lam if lam else float('nan'):.4f})")
    return lam, psi, pdata

print("="*100)
run(30, 8)
print("="*100)
run(13, 8)
print("="*100)
run(30, 4)
