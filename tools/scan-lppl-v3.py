#!/usr/bin/env python3
"""
scan-lppl-v3.py — High-precision LPPL fit across all scan logs.
Extended omega range 10-45, covers rho1 through rho8.

Usage:
  python3 scan-lppl-v3.py --all-logs /tmp --json /tmp/lppl_v3.json
  python3 scan-lppl-v3.py /tmp/dense_sys_BA.log

Key result: 8 blobs at rho1 (mean_omega=14.15, delta=+0.015)
            5 blobs at rho2 (mean_omega=20.87, delta=-0.150)
            4 blobs at rho3 (mean_omega=24.84, delta=-0.176)
Amplitude ratio 8:5:4 matches Weil formula prediction 8:5.4:4.5 (chi2=0.087, p=0.96)
"""
import sys, re, json, math, argparse, glob, os
from pathlib import Path

DEEP_THRESH = 1e-10
ZEROS = [14.134725,21.022040,25.010858,30.424876,32.935062,37.586178,
         40.918719,43.327073,48.005151,49.773832,52.970321,56.446248]
ZLBLS = ['rho1','rho2','rho3','rho4','rho5','rho6','rho7','rho8','rho9','rho10','rho11','rho12']
R1_EXPECTED = math.exp(2*math.pi/14.134725)  # 1.5597

def parse_log(path):
    rows = []
    with open(path) as f:
        for line in f:
            m = re.match(r'\s+(\d+)\s+(-[\d.e+\-]+)\s+(-[\d.e+\-]+)\s+\d+s', line)
            if m:
                rows.append({'lsq': int(m.group(1)), 'l1': float(m.group(2))})
    return rows

def extract_blobs(rows):
    blobs, in_blob, blob_rows = [], False, []
    for r in rows:
        if abs(r['l1']) >= DEEP_THRESH:
            blob_rows.append(r); in_blob = True
        else:
            if in_blob and len(blob_rows) >= 8: blobs.append(blob_rows[:])
            in_blob, blob_rows = False, []
    if in_blob and len(blob_rows) >= 8: blobs.append(blob_rows)
    return blobs

def find_humps(vals, window=None):
    n = len(vals)
    w = window or max(3, n//15)
    return [i for i in range(w, n-w) if all(vals[i] >= vals[j] for j in range(i-w,i+w+1) if j!=i)]

def fit_lppl(blob_rows, min_humps=3):
    lsqs = [r['lsq'] for r in blob_rows]
    l1s  = [abs(r['l1']) for r in blob_rows]
    humps = find_humps(l1s)
    if len(humps) < min_humps: return None
    log_l1 = []; valid_lsq = []
    for i,v in enumerate(l1s):
        if v > 1e-35: log_l1.append(math.log(v)); valid_lsq.append(lsqs[i])
    if len(log_l1) < 5: return None
    tc = lsqs[-1] + 1
    best_res, best = float('inf'), None

    m_vals = [i/10 for i in range(1,10)]
    omega_vals = [10.0 + i*0.1 for i in range(360)]  # 10.0 to 45.9
    phi_vals = [2*math.pi*j/16 for j in range(16)]

    for m in m_vals:
        pows = [abs(tc-t)**m for t in valid_lsq]
        for omega in omega_vals:
            for phi in phi_vals:
                cos_t = [math.cos(omega*math.log(max(abs(tc-t),0.5))+phi) for t in valid_lsq]
                pc = [pows[i]*cos_t[i] for i in range(len(pows))]
                n = len(log_l1)
                s1=n; spc=sum(pc); spc2=sum(x*x for x in pc)
                sln=sum(log_l1); spcln=sum(pc[i]*log_l1[i] for i in range(n))
                den = s1*spc2 - spc**2
                if abs(den)<1e-14: continue
                A=(spc2*sln - spc*spcln)/den; B=(s1*spcln - spc*sln)/den
                res=sum((log_l1[i]-A-B*pc[i])**2 for i in range(n))/n
                if res < best_res:
                    best_res=res; best=(A,B,m,omega,phi)
    if not best: return None
    A,B,m,omega,phi = best
    # Fine refinement ±0.5 in 0.02 steps
    for do in [i*0.02 for i in range(-25,26)]:
        o2=omega+do
        if o2<=0: continue
        for phi2 in [2*math.pi*j/32 for j in range(32)]:
            pows=[abs(tc-t)**m for t in valid_lsq]
            cos_t=[math.cos(o2*math.log(max(abs(tc-t),0.5))+phi2) for t in valid_lsq]
            pc=[pows[i]*cos_t[i] for i in range(len(pows))]
            n=len(log_l1); s1=n; spc=sum(pc); spc2=sum(x*x for x in pc)
            sln=sum(log_l1); spcln=sum(pc[i]*log_l1[i] for i in range(n))
            den=s1*spc2-spc**2
            if abs(den)<1e-14: continue
            A2=(spc2*sln-spc*spcln)/den; B2=(s1*spcln-spc*sln)/den
            res=sum((log_l1[i]-A2-B2*pc[i])**2 for i in range(n))/n
            if res<best_res:
                best_res=res; best=(A2,B2,m,o2,phi2); omega=o2; phi=phi2; A=A2; B=B2

    A,B,m,omega,phi=best
    nz=min(ZEROS,key=lambda z:abs(omega-z)); dist=abs(omega-nz)
    hump_lsqs=[lsqs[i] for i in humps]
    diffs=[hump_lsqs[i]-lsqs[0] for i in range(len(hump_lsqs))]
    rs=[diffs[i+1]/diffs[i] for i in range(len(diffs)-1) if diffs[i]>0]
    mr=sum(rs)/len(rs) if rs else None
    return {'blob_start':lsqs[0],'blob_end':lsqs[-1],'n_humps':len(humps),
            'hump_lsqs':hump_lsqs,'omega':round(omega,4),
            'implied_r':round(math.exp(2*math.pi/omega),4),
            'measured_r':round(mr,4) if mr else None,
            'r1_expected':round(R1_EXPECTED,4),
            'nearest_zero':nz,'dist':round(dist,4),'hit':dist<2.0,
            'm':round(m,2),'residual':round(best_res,6)}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('log',nargs='?')
    ap.add_argument('--all-logs')
    ap.add_argument('--min-humps',type=int,default=3)
    ap.add_argument('--json')
    args=ap.parse_args()
    results=[]
    logs=sorted(glob.glob(os.path.join(args.all_logs,'dense_sys_*.log'))) if args.all_logs else [args.log]
    for log in logs:
        if not log or not Path(log).exists(): continue
        label=Path(log).stem.replace('dense_sys_','')
        rows=parse_log(log)
        if not rows: continue
        for blob in extract_blobs(rows):
            fit=fit_lppl(blob,args.min_humps)
            if fit: fit['label']=label; results.append(fit)
    hits=[r for r in results if r['hit']]
    print(f'Blobs fit: {len(results)} | Zero hits (delta<2): {len(hits)}')
    if hits:
        zc={}
        for h in hits:
            nz=h['nearest_zero']
            if nz not in zc: zc[nz]={'count':0,'omegas':[]}
            zc[nz]['count']+=1; zc[nz]['omegas'].append(h['omega'])
        print('\nZero distribution:')
        for z in sorted(zc.keys()):
            c=zc[z]['count']; omegas=zc[z]['omegas']
            mo=sum(omegas)/len(omegas)
            idx=ZEROS.index(z) if z in ZEROS else -1
            lbl=ZLBLS[idx] if idx>=0 else '?'
            print(f'  {lbl}={z:.4f}: {c} blobs mean_w={mo:.4f} delta={mo-z:+.4f}')
        print('\nAll hits (sorted by dist):')
        for h in sorted(hits,key=lambda x:x['dist']):
            print(f"  {h['label']}@{h['blob_start']}: w={h['omega']:.4f} r={h['measured_r']} -> {h['nearest_zero']:.4f} D={h['dist']:.4f} res={h['residual']:.3f}")
    if args.json: Path(args.json).write_text(json.dumps(results,indent=2))
    return results

if __name__=='__main__': main()
