# D(t) MACHINE (Cor 4.1, Suzuki [14]).  Two validations + the defect derivation made concrete:
# (A) g_xi(t) from EQ_403 (prime/Gamma/Lerch, UNCONDITIONAL) == spectral -2*sum_{g>0}(1-cos tg)/g^2 (needs zeros).
# (B) illustrate D_rho(t): a HYPOTHETICAL off-line zero forces an e^{2*delta*t} excess in (1/2pi)||S_t||^2.
import mpmath as mp
mp.mp.dps=30
PSI14_MINUS_LOGPI = mp.digamma(mp.mpf(1)/4) - mp.log(mp.pi)     # = Omega(0) = -5.372...
TRIGAMMA14 = mp.zeta(2, mp.mpf(1)/4)                            # Phi(1,2,1/4)=sum 1/(n+1/4)^2 = psi'(1/4)
def vonmangoldt_upto(N):
    Lam=[mp.mpf(0)]*(N+1); is_c=[False]*(N+1)
    for p in range(2,N+1):
        if not is_c[p]:
            lp=mp.log(p); pk=p
            while pk<=N:
                Lam[pk]=lp; pk*=p
            for m in range(2*p,N+1,p): is_c[m]=True
    return Lam
def g_xi(t, Lam):
    t=mp.mpf(t); et=mp.e**t
    s = -4*(mp.e**(t/2)+mp.e**(-t/2)-2)
    N=int(mp.floor(et))
    for n in range(2,N+1):
        if Lam[n]>0: s += Lam[n]/mp.sqrt(n)*(t-mp.log(n))
    s += -(t/2)*PSI14_MINUS_LOGPI
    s += -mp.mpf(1)/4*(TRIGAMMA14 - mp.e**(-t/2)*mp.lerchphi(mp.e**(-2*t),2,mp.mpf(1)/4))
    return s
def g_spectral(t, gammas):
    t=mp.mpf(t)
    return -2*mp.fsum((1-mp.cos(t*g))/g**2 for g in gammas)
# load zeros
G=[]
with open('data/zeros_100k.txt') as f:
    for line in f:
        line=line.strip()
        if line:
            try: G.append(mp.mpf(line))
            except: pass
        if len(G)>=5000: break
print(f"loaded {len(G)} zeros; g_xi[EQ_403 prime form] vs g_spectral[zero form]:")
Lam=vonmangoldt_upto(200000)
print(f"{'t':>4} {'g_xi (primes)':>18} {'g_spectral (zeros)':>20} {'diff':>12}")
for t in [0.5,1.0,2.0,3.0,4.0]:
    a=g_xi(t,Lam); b=g_spectral(t,G)
    print(f"{t:>4} {mp.nstr(a,10):>18} {mp.nstr(b,10):>20} {mp.nstr(a-b,3):>12}")
# (B) off-line excess: mass_gamma(t)=(e^{-i gamma t}-1)/gamma; off-line zero => gamma0 = -Im + i*delta, delta=beta-1/2
print("\n(B) off-line-zero EXCESS |mass(t)|^2 = |e^{-i g0 t}-1|^2/|g0|^2, g0=-100 + i*delta (delta=beta-1/2):")
print(f"{'t':>4} " + "".join(f"delta={d}".rjust(14) for d in [0.0,0.1,0.2,0.3]))
for t in [1,3,5,8]:
    row=f"{t:>4} "
    for d in [0.0,0.1,0.2,0.3]:
        g0=mp.mpf(-100)+1j*mp.mpf(str(d))
        m=abs((mp.e**(-1j*g0*t)-1)/g0)**2
        row+=mp.nstr(m,4).rjust(14)
    print(row)
print("delta=0 (on-line): bounded oscillation; delta>0: grows ~ e^{2*delta*t} => D(t)->inf, cannot be identically 0.")
