#!/usr/bin/env python3
"""
scan-lppl-convergence.py — Convergence analysis of LPPL omega measurements.

Reads LPPL fit results and shows:
1. Running mean of omega vs each Riemann zero as blobs accumulate
2. Amplitude ordering test vs Weil formula prediction
3. Chi-squared fit
4. Hump-ratio test

Usage: python3 scan-lppl-convergence.py /tmp/lppl_v3_latest.json
"""
import json, math, sys
from pathlib import Path

ZEROS = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178, 40.918719, 43.327073]
ZLBLS = ['rho1','rho2','rho3','rho4','rho5','rho6','rho7','rho8']
R1_THEORY = math.exp(2*math.pi/14.134725)

def get_cluster(data, z, delta=2.0, max_res=20, min_omega=11):
    return sorted([d for d in data
                   if abs(d.get('omega',0)-z) < delta
                   and d.get('residual',999) < max_res
                   and d.get('omega',0) > min_omega],
                  key=lambda x: x['blob_start'])

def main():
    f = sys.argv[1] if len(sys.argv)>1 else '/tmp/lppl_v3_latest.json'
    data = json.loads(Path(f).read_text())
    print(f"Total blobs: {len(data)}\n")

    # Convergence for each detected zero
    for z, lbl in zip(ZEROS[:4], ZLBLS[:4]):
        cluster = get_cluster(data, z)
        if not cluster:
            print(f"{lbl}={z:.4f}: no hits\n")
            continue
        n = len(cluster)
        omegas = [d['omega'] for d in cluster]
        mean = sum(omegas)/n
        std = math.sqrt(sum((o-mean)**2 for o in omegas)/(n-1)) if n>1 else 0
        se = std/math.sqrt(n)
        ci = (mean-1.96*se, mean+1.96*se)
        t = (mean-z)/se if se>0 else 0
        print(f"{lbl}={z:.6f}: n={n} mean_w={mean:.6f} std={std:.4f} SE={se:.4f}")
        print(f"  delta={mean-z:+.6f} ({(mean-z)/z*100:+.4f}%)")
        print(f"  t-stat={t:.4f}  95%CI=[{ci[0]:.4f},{ci[1]:.4f}]  rho_inside={ci[0]<=z<=ci[1]}")
        print(f"  Blobs: lsq {cluster[0]['blob_start']} to {cluster[-1]['blob_start']}")
        print(f"  Running convergence:")
        running = []
        for i, d in enumerate(cluster):
            running.append(d['omega'])
            m = sum(running)/len(running)
            print(f"    N={i+1:>2}: lsq={d['blob_start']:>7} w={d['omega']:>7.3f} mean={m:>9.5f} delta={m-z:>+9.6f}")
        print()

    # Amplitude ordering chi-squared test
    print("Weil amplitude fit:")
    zc = {}
    quality = [d for d in data if d.get('omega',0)>11 and d.get('residual',999)<20 and d.get('hit')]
    for d in quality: zc[d['nearest_zero']]=zc.get(d['nearest_zero'],0)+1
    n1 = zc.get(ZEROS[0], 0)
    if n1==0: print("No rho1 hits found."); return
    chi2=0; dof=0
    print(f"  (n_rho1={n1})")
    print(f"  {'rho':>7} {'obs':>5} {'pred':>7} {'O/P':>6} {'chi2':>8}")
    detected_zeros = [z for z in ZEROS[:6] if zc.get(z,0)>0]
    for z in ZEROS[:6]:
        o = zc.get(z,0)
        p = n1 * ZEROS[0]/z
        ct = (o-p)**2/p
        chi2 += ct; dof += 1
        print(f"  {z:>7.3f} {o:>5} {p:>7.2f} {o/p if p>0 else 0:>6.3f} {ct:>8.4f}")
    print(f"  chi2={chi2:.4f} ({dof}dof)  p>0.05: {chi2<12.59}")

    # Hump ratio test for rho1 hits
    print("\nHump-ratio test (theory r1=1.5597):")
    rho1 = [d for d in quality if d['nearest_zero']==ZEROS[0]]
    for d in rho1:
        mr = d.get('measured_r')
        ir = math.exp(2*math.pi/d['omega'])
        if mr:
            print(f"  {d['label']}@{d['blob_start']}: w={d['omega']} measured_r={mr:.4f} implied={ir:.4f} theory={R1_THEORY:.4f} delta={abs(mr-R1_THEORY):.4f}")

if __name__=='__main__': main()
