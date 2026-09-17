import mpmath as mp
mp.mp.dps = 30
G = [mp.mpf(x) for x in open("operators/_zeros150.txt").read().split()[:12]]

def eta(s): return (1-2**(1-s))*mp.zeta(s)   # Psi_s(0) = -f*Gamma(s)*eta(s); Gamma!=0 => BC zero <=> eta(s)=0

print("=== (B) roots of Psi_s(0)=0  [i.e. eta(s)=0], seeded at 1/2+i*gamma_k, in FULL complex plane ===")
print("    If the construction is faithful, each root has Re(s)=1/2 exactly and Im=gamma_k.")
maxdev = mp.mpf(0)
for gk in G:
    root = mp.findroot(eta, mp.mpc('0.5', gk))
    dev = abs(root.real - mp.mpf('0.5'))
    maxdev = max(maxdev, dev)
    print(f"  seed 1/2+{float(gk):.5f}i -> root Re={float(root.real):.12f}  Im={float(root.imag):.6f}  |eta(root)|={float(abs(eta(root))):.1e}")
print(f"  --> max |Re(root)-1/2| over 12 zeros: {float(maxdev):.2e}   (eigenvalue i(1/2-s) => real gamma)")

print("\n=== (C) spurious eta-zeros from (1-2^{1-s})=0 : s=1+2*pi*i*k/ln2  (eta=0, zeta!=0) ===")
for k in [1,2,3]:
    s = mp.mpc(1, 2*mp.pi*k/mp.log(2))
    print(f"  k={k}: s={complex(s):.4f}  |1-2^(1-s)|={float(abs(1-2**(1-s))):.2e}  |zeta(s)|={float(abs(mp.zeta(s))):.4f}  |eta(s)|={float(abs(eta(s))):.1e}")
print("  => these are eigenvalue labels of R_hat that are NOT zeta zeros (Re=1, off critical line):")
print("     the construction's known spurious branch. Yakaboylu excises them via the domain D(R_hat_c).")
