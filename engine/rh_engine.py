#!/usr/bin/env python3
"""THE RH ENGINE — unstoppable, learning, creative solving loop.
Invariant: the ONLY halt is PROOF_FOUND. Low queue -> generate. Stagnation -> escalate.
Elimination -> fuel. Never "exhausted." See ENGINE.md."""
import json, os, subprocess, time, sys, datetime, hashlib
# HARD GUARD (Jacob 2026-07-18): NEVER bill the Anthropic API key — that is real cash owed.
# Scrub it so no subprocess (claude -p etc.) can ever use it. Max subscription + FREE routes only.
os.environ.pop("ANTHROPIC_API_KEY",None); os.environ.pop("ANTHROPIC_AUTH_TOKEN",None)

# portable ROOT: Mac (Sam) or Olympus (the always-on daemon home)
ROOT="/Users/scones/projects/riemann" if os.path.exists("/Users/scones") else os.path.expanduser("~/projects/riemann")
ENG=f"{ROOT}/engine"

# LIT-GATE WIRING (2026-07-28): the heuristic keyword judge below (classify()) has no idea
# whether a "LEAD" is actually already-published territory or a known-dead route — that lit-gate
# check existed only as arsenal/t01_route_prefilter, a separate tool nobody was forced to run
# against engine-generated methods. Wire it in here so a false LEAD from keyword coincidence
# can never bump bandit yield on rediscovered/dead ground; see arsenal/t12_roundtable_gate.py
# for the sibling fix on the human-run Roundtable side.
sys.path.insert(0, f"{ROOT}/arsenal")
try:
    from t01_route_prefilter import _prior_art_hits
except Exception:
    _prior_art_hits = None
REG=f"{ENG}/methods.json"
LOG=f"{ENG}/outcomes.jsonl"
PROOF=f"{ROOT}/results/PROOF_FOUND.marker"
def _py(v):  # use the venv if present (Mac), else system python3 (Olympus has numpy/mpmath/scipy)
    p=f"{ROOT}/{v}_venv/bin/python"; return p if os.path.exists(p) else "python3"
PYS={"cvs":_py("cvs"),"spectral":_py("spectral"),"flint":_py("flint"),"sys":"python3"}
MIN_QUEUE=4            # refill trigger
MIN_GEN_INTERVAL=600   # 10min between generations — FREE via Olympus Ollama, so keep the GPU productive
                       # (paid Claude is capped SEPARATELY in budget.py regardless of this)
STAGNATION=5          # consecutive nulls -> escalate
GEN_BATCH=4           # methods per generation

def now(): return datetime.datetime.now().isoformat(timespec="seconds")
def ntfy_breakthrough(msg):
    """Phone ping ONLY for a gate-passing breakthrough. Nothing else reaches the phone."""
    try:
        subprocess.run(["node",os.path.expanduser("~/tools/notify/notify.js"),msg.upper(),
            "--title","RH BREAKTHROUGH","--priority","max"],
            capture_output=True,timeout=25)
    except Exception: pass
def load(p,d): 
    try: return json.load(open(p))
    except: return d
def save(p,o): json.dump(o,open(p,"w"),indent=2)
def logline(o): open(LOG,"a").write(json.dumps(o)+"\n")

def seed_registry():
    """Seed from this session's real runnable scripts + known dead classes."""
    if os.path.exists(REG): return load(REG,{})
    m={"methods":[],"eliminated_classes":[
        "local/factorization positivity (moment/operator/SUSY/TP all show positivity is product-emergent, A14)",
        "linear-algebra reshuffles of tau (pivots/eigenvector/Schur/H-matrix/kernel-count, A1/A2/A9/A10/A11)",
        "zero-statistics (GUE/moments/gaps/FHK: all confirm, none crack — forcing is analytic not statistical)",
        "pointwise majorant of prime side (Beurling-Selberg), large sieve, de Branges naive (OBSTRUCTION_MAP)"],
      "families":{}, "consecutive_nulls":0, "cycles":0}
    seeds=[
      ("A6_cross_L","cross_L_A6.py","universality","done"),
      ("A12_susy","susy_factorization_probe.py","susy","lead"),
      ("A13_susy_basis","susy_basis_A13.py","susy","lead"),
      ("A14_toeplitz_tp","toeplitz_tp_test.py","total_positivity","done"),
    ]
    for name,scr,fam,st in seeds:
        m["methods"].append({"name":name,"code":f"operators/{scr}","venv":"cvs","family":fam,
                             "status":st,"yield":1.0 if st=="lead" else 0.3,"born":now()})
    save(REG,m); return m

def classify(out, err, rc):
    """JUDGE (heuristic core; LLM judge is the escalation). PROOF/LEAD/NULL/WALL."""
    t=(out or "")+(err or "")
    if os.path.exists(PROOF) or "PROOF_FOUND" in t: return "PROOF"
    if rc!=0 and "Traceback" in t: return "ERROR"
    low=t.lower()
    lead_sig=["strictly positive","new identity","closed form","converges to","factor","clean","smooth",
              "gap","structure","0 errors","certified positive"]
    wall_sig=["nonlocal","not manifest","fails","no clean","no anomaly","consistent with","degenerate",
              "not an h-matrix","not converging","emergent"]
    ls=sum(s in low for s in lead_sig); ws=sum(s in low for s in wall_sig)
    if ls>ws and ls>=2: return "LEAD"
    if ws>=2: return "WALL"
    return "NULL"

def prior_art_override(mth):
    """Cross-check name/family/rationale/novelty text against the arsenal prior-art DB.
    Returns (override_class, hit) or (None, None). A DEAD/WALL hit forces WALL regardless of
    the heuristic judge's keyword count (stops re-crediting a known-dead route as a LEAD).
    A KNOWN hit (published, non-dead — e.g. Jensen/GORTTW territory) forces KNOWN_TERRITORY so
    rediscovering a published result never bumps bandit yield as if it were fresh structure."""
    if _prior_art_hits is None:
        return None, None
    text=" ".join([mth.get("name",""),mth.get("family",""),
                    mth.get("rationale",""),mth.get("novelty","")])
    try:
        hits=_prior_art_hits("",text)
    except Exception:
        return None, None
    dead_wall=[h for h in hits if h[1] in("DEAD","WALL")]
    known=[h for h in hits if h[1]=="KNOWN"]
    if dead_wall: return "WALL", dead_wall[0]
    if known: return "KNOWN_TERRITORY", known[0]
    return None, None

def run_method(mth,reg):
    code=f"{ROOT}/{mth['code']}"
    if not os.path.exists(code):
        return {"status":"NO_CODE","note":"code missing — flag for generator to build"}
    py=PYS.get(mth.get("venv","sys"),PYS["sys"])
    print(f"[{now()}] RUN {mth['name']} ({mth['family']}) via {mth['venv']}",flush=True)
    try:
        r=subprocess.run([py,code],capture_output=True,text=True,timeout=1200,cwd=ROOT)
        cls=classify(r.stdout,r.stderr,r.returncode)
    except subprocess.TimeoutExpired:
        r=None; cls="TIMEOUT"
    tail=(r.stdout[-1500:] if r else "")+("\n[stderr]"+r.stderr[-400:] if r and r.stderr else "")
    out={"ts":now(),"method":mth["name"],"family":mth["family"],"class":cls,"tail":tail}
    logline(out)
    print(f"[{now()}]   -> {cls}",flush=True)
    return out

def generate(reg):
    """CREATIVITY CORE. Calls Claude with full state -> new methods (+code). Never runs dry.
    THROTTLED so the loop runs forever on free compute; Claude generation only every MIN_GEN_INTERVAL."""
    last=reg.get("last_gen_ts",0)
    if time.time()-last < MIN_GEN_INTERVAL:
        print(f"[{now()}] generation throttled ({int((MIN_GEN_INTERVAL-(time.time()-last))/60)}m to next); loop continues",flush=True); return
    reg["last_gen_ts"]=time.time()
    print(f"[{now()}] GENERATE (queue low / stagnation) — invoking creativity core",flush=True)
    try:
        r=subprocess.run(["python3",f"{ENG}/generate_methods.py"],capture_output=True,text=True,timeout=900)
        print(r.stdout[-800:],flush=True)
    except Exception as e:
        print(f"[{now()}] generator unavailable ({e}); wrote prompt for next session",flush=True)

def cycle():
    reg=seed_registry()
    reg["cycles"]=reg.get("cycles",0)+1
    if os.path.exists(PROOF):
        print(f"[{now()}] *** PROOF_FOUND — machine halts (the one allowed exit). ***");save(REG,reg);return "PROOF"
    pending=[m for m in reg["methods"] if m["status"] in("pending","lead") and os.path.exists(f"{ROOT}/{m['code']}")]
    pending.sort(key=lambda m:-m.get("yield",0.5))
    if not pending:
        generate(reg); save(REG,reg); return "GENERATED"
    mth=pending[0]
    res=run_method(mth,reg)
    # LIT-GATE: force-correct the heuristic judge against known prior art before LEARN.
    pa_class,pa_hit=prior_art_override(mth)
    if pa_class and res.get("class") in("LEAD","PROOF"):
        logline({"ts":now(),"method":mth["name"],"lit_gate_override":
                 f"{res['class']} -> {pa_class} (prior-art: {pa_hit[0]} [{pa_hit[1]}] {pa_hit[2][:100]})"})
        print(f"[{now()}]   LIT-GATE OVERRIDE: {res['class']} -> {pa_class} (matched '{pa_hit[0]}')",flush=True)
        res["class"]=pa_class
    # LEARN
    if res.get("class")=="PROOF":
        # VERIFICATION GATE: a claim is REPORT until an independent rigorous check makes it HANDLE.
        try:
            import verify_gate; handle,report=verify_gate.verify()
        except Exception as e:
            handle,report=False,f"gate error: {e}"
        logline({"ts":now(),"method":mth["name"],"gate":report,"handle":handle})
        if handle:
            open(PROOF,"a").write(f"{mth['name']} {now()} | {report}\n"); mth["status"]="PROOF"
            ntfy_breakthrough(f"VERIFIED proof-candidate from {mth['name']}!! independently confirmed. CHECK NOW.")
        else:
            mth["status"]="unverified_claim"   # downgraded — NO halt, NO phone ping, keep running
            print(f"[{now()}]   PROOF claim REJECTED by gate: {report}",flush=True)
    elif res.get("class")=="LEAD":
        mth["status"]="explored"; mth["yield"]=mth.get("yield",0.5)+0.5; reg["consecutive_nulls"]=0
    elif res.get("class")=="WALL":
        mth["status"]="wall"; reg.setdefault("eliminated_classes",[]).append(f"{mth['family']} (via {mth['name']}, {now()})")
        reg["consecutive_nulls"]=reg.get("consecutive_nulls",0)+1
    else:
        mth["status"]="explored"; reg["consecutive_nulls"]=reg.get("consecutive_nulls",0)+1
    fam=reg.setdefault("families",{}); fam[mth["family"]]=fam.get(mth["family"],0)+1
    # anti-stall: refill / escalate
    live=[m for m in reg["methods"] if m["status"] in("pending","lead")]
    if len(live)<MIN_QUEUE or reg.get("consecutive_nulls",0)>=STAGNATION:
        generate(reg)
        if reg.get("consecutive_nulls",0)>=STAGNATION:
            print(f"[{now()}] STAGNATION={reg['consecutive_nulls']} -> ESCALATE (force different method KIND)",flush=True)
            reg["consecutive_nulls"]=0
    save(REG,reg); return res.get("class","NULL")

if __name__=="__main__":
    n=int(sys.argv[1]) if len(sys.argv)>1 else 1
    for _ in range(n):
        if cycle()=="PROOF": break
