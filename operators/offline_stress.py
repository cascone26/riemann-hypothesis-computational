"""Off-line stress test of screw function g ([13] converse mechanism). Vectorized."""
import numpy as np
G = np.load('zeros2000.npy').astype(float)

def g_vec(zset, ts, chunk=400):
    zr = zset.real[None,:]; zi = zset.imag[None,:]; z2 = (zset*zset)[None,:]
    out = np.empty(len(ts))
    for i in range(0, len(ts), chunk):
        t = ts[i:i+chunk][:,None]
        a = np.exp(1j*zset[None,:]*t); b = np.exp(-1j*zset[None,:]*t)
        out[i:i+chunk] = -0.5*np.real(np.sum((a-1)*(b-1)/z2, axis=1))
    return out

def build_set(g0, eps):
    base = np.concatenate([G, -G]).astype(complex)
    if eps == 0.0: return base
    base = base[np.abs(np.abs(base.real)-g0) > 1e-9]
    quad = np.array([g0+1j*eps, g0-1j*eps, -g0+1j*eps, -g0-1j*eps])
    return np.concatenate([base, quad])

if __name__ == "__main__":
    idx=100; g0=G[idx]
    print(f"gamma0={g0:.4f}; g(t)>0 somewhere <=> RH-false signature. Predicted crossing at eps*t~7.")
    print(f"{'eps':>6} {'first t: g>0':>14} {'sup g':>12} {'t@sup':>9} {'eps*t_cross':>11}")
    for eps in [0.05,0.1,0.2,0.4]:
        Tmax = min(400.0, 12.0/eps + 60)   # enough to see crossing
        ts = np.linspace(0.5, Tmax, int(min(400000, Tmax*2500)))
        gv = g_vec(build_set(g0,eps), ts)
        pos = np.where(gv>0)[0]
        supg=gv.max(); tsup=ts[np.argmax(gv)]
        if len(pos):
            tf=ts[pos[0]]; print(f"{eps:6.2f} {tf:14.3f} {supg:12.4e} {tsup:9.3f} {eps*tf:11.3f}")
        else:
            print(f"{eps:6.2f} {'(none<=%.0f)'%Tmax:>14} {supg:12.4e} {tsup:9.3f} {'--':>11}")
