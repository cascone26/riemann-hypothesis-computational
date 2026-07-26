#!/usr/bin/env python3
"""S1 Onset Lemma probe: how does a prime power's contribution P_q(k,l) to the even Weil block scale as
eps = a - (1/2)log q -> 0+ (its onset), and is P_q RANK-ONE (a manifest signed square)?

If P_q ~ eps^n * (u u^T) then: onset order = n, and the drain is a manifest square times an arithmetic
weight Lambda(q)/sqrt(q) -> the 'race' has an estimate-IMMUNE (SOS) drain term. This is the keystone of the
per-gap race reformulation. Uses the II/dterms float machinery verified vs Arb in relay_decomp (§218).
"""
import numpy as np, cmath, math

def Pp(p,t,s):
    if p==0: return 1/s
    if p==1: return t/s - 1/s**2
    return t*t/s - 2*t/s**2 + 2/s**3
def II(p,alpha,omega,sc,t0,t1):
    s=complex(-alpha,omega)
    def F(tt):
        t=complex(tt); return cmath.exp(s*t)*Pp(p,t,s)
    v=F(t1)-F(t0)
    return v.real if sc=='cos' else v.imag
def dterms(k,l,A,B):
    if k==l:
        return [(2*A*k*k*B*B,0,k*B,'cos'),(-k*k*B*B,1,k*B,'cos'),(-k*B,0,k*B,'sin')]
    al=k*l*B*(1.0/(k-l)-1.0/(k+l)); ak=-k*l*B*(1.0/(k-l)+1.0/(k+l))
    return [(al,0,l*B,'sin'),(ak,0,k*B,'sin')]

def P_q(q, a, modes):
    """the q-prime-power contribution to the even block, at scale a."""
    A=a; B=math.pi/(2*A); TWOA=2*A; LOGQ=math.log(q); Wq=math.log(_rad(q))/math.sqrt(q)
    def entry(k,l):
        tot=0.0
        for (ad,pd,om,sc) in dterms(k,l,A,B):
            tot+=Wq*ad*II(pd+1,0.0,om,sc,LOGQ,TWOA)+(-Wq*LOGQ)*ad*II(pd,0.0,om,sc,LOGQ,TWOA)
        return tot
    M=np.array([[entry(k,l) for l in modes] for k in modes])
    return (M+M.T)/2
def _rad(q):  # radical: the prime p for prime power q=p^m
    p=2
    while q% p: p+=1
    return p

modes=[1,3,5,7]
print(f"S1 ONSET LEMMA — universality across prime powers, even modes {modes}")
print("Claim: P_q(a) ~ -Lambda(q)/sqrt(q) * kappa * eps^3 * (l x l),  l_k = k,  eps=a-½log q")
print(f"\n{'q':>4} {'onset a0':>10} {'order n':>8} {'|lam2|/|lam1|':>13} {'u_k/k (should be const)':>34} {'sign':>6}")
for q in [3,4,5,7,9]:
    a0=0.5*math.log(q)
    data=[]
    for eps in [0.0005,0.001,0.002,0.004,0.008]:
        M=P_q(q,a0+eps,modes); data.append((eps,np.linalg.norm(M)))
    eps=np.array([d[0] for d in data]); nf=np.array([d[1] for d in data])
    n,_=np.polyfit(np.log(eps),np.log(nf),1)
    M=P_q(q,a0+0.004,modes); w,V=np.linalg.eigh(M); idx=np.argsort(np.abs(w))[::-1]
    r=abs(w[idx[1]])/abs(w[idx[0]]); u=V[:,idx[0]]
    if u[0]<0: u=-u  # sign convention
    uk_over_k=u/np.array(modes,float)
    sign='NEG' if w[idx[0]]<0 else 'POS'
    print(f"{q:>4} {a0:>10.5f} {n:>8.3f} {r:>13.2e} {np.array2string(uk_over_k,precision=4):>34} {sign:>6}")

# extract kappa (universal constant) from q=3
q=3; a0=0.5*math.log(q); eps=0.004
M=P_q(q,a0+eps,modes); w=np.linalg.eigvalsh(M); lam1=w[np.argmax(np.abs(w))]
Wq=math.log(3)/math.sqrt(3)
# lam1 = -Wq*kappa*eps^3*|l|^2 ; |l|^2 = sum k^2 for k in modes
l2=sum(k*k for k in modes)
kappa=-lam1/(Wq*eps**3*l2)
print(f"\nq=3: lam1={lam1:.6e}, Wq=Lambda(3)/sqrt(3)={Wq:.5f}, |l|^2=sum k^2={l2}")
print(f"  => kappa = -lam1/(Wq*eps^3*|l|^2) = {kappa:.6f}   (the universal onset constant)")
print("READ: if order n≈3, |lam2/lam1|≈0, u_k/k constant, sign NEG for ALL q => ONSET LEMMA holds universally:")
print("  each prime power drains as a rank-one NEGATIVE square in the FIXED direction l_k=k, weight Lambda(q)/sqrt(q).")
