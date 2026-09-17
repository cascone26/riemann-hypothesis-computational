import mpmath as mp
mp.mp.dps = 25
# Verify DH construction: D(s) = (1-iκ)/2 L(s,χ) + (1+iκ)/2 L(s,χ̄), χ mod 5, χ(2)=i.
# L(s,χ) = 5^-s Σ_{r=1}^5 χ(r) ζ(s, r/5) (Hurwitz). Find zeros; confirm OFF-LINE ones exist (~t=85).
kap = (mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
chi = {1:1, 2:1j, 3:-1j, 4:-1, 0:0}   # χ(n)=χ(n mod 5)
def L(s, conj=False):
    tot = mp.mpf(0)
    for r in range(1,6):
        c = chi[r%5]
        if conj: c = mp.conj(c)
        tot += c * mp.zeta(s, mp.mpf(r)/5)
    return 5**(-s) * tot
def D(s):
    return (1-1j*kap)/2 * L(s,False) + (1+1j*kap)/2 * L(s,True)

print(f"κ={mp.nstr(kap,8)}")
# sanity: D should satisfy functional eq via ξ(s)=(π/5)^{-s/2}Γ((1+s)/2)D(s), ξ(s)=ξ(1-s)
def xi(s): return (mp.pi/5)**(-s/2)*mp.gamma((1+s)/2)*D(s)
s0=mp.mpf('0.3')+17j
print(f"functional eq check |ξ(s)-ξ(1-s)|/|ξ(s)| at s=0.3+17i: {mp.nstr(abs(xi(s0)-xi(1-s0))/abs(xi(s0)),4)}  (should be ~0)")

# find zeros on critical line (Re=1/2) low, and search for OFF-LINE zeros near t~85
print("\nzeros of D on Re=1/2 (first few, sign changes of Re/Im of D(1/2+it)):")
prev=D(mp.mpf('0.5')+10j); onl=[]
t=mp.mpf('10')
while t<30 and len(onl)<4:
    v=D(mp.mpf('0.5')+1j*t)
    if mp.re(prev)*mp.re(v)<0 or mp.im(prev)*mp.im(v)<0:
        try:
            z=mp.findroot(lambda x: D(mp.mpf('0.5')+1j*x), float(t))
            onl.append(float(z))
        except: pass
    prev=v; t+=mp.mpf('0.5')
print(f"  on-line t: {[round(z,3) for z in onl]}")

# OFF-LINE search: |D(σ+it)| minima for σ≠1/2 near t=80..88 (literature: DH off-line zero ~ 0.808+85.7i region)
print("\noff-line zero search (|D| minima off Re=1/2, t in 80..90):")
found=[]
for sig in [mp.mpf('0.7'),mp.mpf('0.8'),mp.mpf('0.9'),mp.mpf('1.0'),mp.mpf('1.1')]:
    for t in [mp.mpf(x)/10 for x in range(800,905,3)]:
        v=abs(D(sig+1j*t))
        if v<0.02:
            found.append((float(sig),float(t),float(v)))
if found:
    found.sort(key=lambda x:x[2])
    for sig,t,v in found[:6]: print(f"  σ={sig} t={t}: |D|={v:.4f}")
else:
    print("  no |D|<0.02 on the coarse grid; refine near candidates")
