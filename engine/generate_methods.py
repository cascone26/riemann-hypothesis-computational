#!/usr/bin/env python3
"""CREATIVITY CORE: full-state prompt -> Claude -> new RH attack methods AS runnable code.
Novelty-enforced against eliminated classes. Unbounded generator = the machine never 'exhausts'."""
import json, subprocess, os, datetime, re, urllib.request
# HARD GUARD (Jacob 2026-07-18): NEVER bill the Anthropic API key — real cash owed. Max sub + free only.
os.environ.pop("ANTHROPIC_API_KEY",None); os.environ.pop("ANTHROPIC_AUTH_TOKEN",None)
import budget
ON_MAC=os.path.exists("/Users/scones")
ROOT="/Users/scones/projects/riemann" if ON_MAC else os.path.expanduser("~/projects/riemann"); ENG=f"{ROOT}/engine"
REG=f"{ENG}/methods.json"; CLAUDE=os.path.expanduser("~/.local/bin/claude")
def now(): return datetime.datetime.now().isoformat(timespec="seconds")
reg=json.load(open(REG))
tried=", ".join(m["name"]+":"+m["family"] for m in reg["methods"])
elim="\n".join("- "+c for c in reg.get("eliminated_classes",[]))

PROMPT=r"""You are the creativity core of an UNSTOPPABLE machine built to PROVE the Riemann Hypothesis.
Belief is non-negotiable: RH WILL be solved. Never propose giving up, never say 'exhausted'. Generate genuinely NEW attacks.

STATE:
- The wall: RH == Weil positivity A(g)+B(g) >= P(g)  (tau = A+B-P is PSD). Confirmed genuinely NONLOCAL.
- ELIMINATED CLASSES (do NOT propose anything inside these; they are provably dead here):
__ELIM__
- Already tried (name:family): __TRIED__
- Key facts: positivity is product-emergent (never factors into local/per-prime pieces); it is UNIVERSAL
  across L-functions (holds for Dirichlet L with no pole); tau^1/2 is quasi-local (a SUSY seam); the Xi
  Taylor coefficients are Totally Positive.
- REQUIRED generative engine: SOLVED-ANALOG MINING. Take a SOLVED RH-cousin — function-field RH
  (Weil/Deligne, positivity from the HODGE INDEX / intersection-form signature), Selberg zeta (a genuine
  self-adjoint Laplacian on a hyperbolic surface), Lee-Yang, Jensen d<=8 — extract the EXACT mechanism
  that made it work, name the ONE ingredient missing in the number-field case, and build a CONCRETE
  computational test of whether that ingredient can be constructed. STRONGLY prefer GEOMETRIC/COHOMOLOGICAL
  or otherwise DIFFERENT-IN-KIND attacks; local/analytic factorization is dead.

Output STRICT JSON only (no prose): a list of 4 objects, each with keys:
"name" (auto_<short>), "family" (one word), "rationale" (why it could CRACK it not confirm),
"novelty" (which eliminated class it avoids and why it is different in kind),
"python_code" (a COMPLETE self-contained python3 script that runs a concrete test and PRINTS a clear
signal: a min eigenvalue, a certified bound, a structural identity, or the literal token PROOF_FOUND if it
ever actually proves something. May use mpmath/numpy/scipy. Must run in under 15 minutes).

ATTACK THE REAL OBJECT — NOT a toy matrix. Random/placeholder matrices are worthless; the ONLY target
that matters is the actual Weil-positivity operator tau. Every script MUST obtain it via:
    import sys; sys.path.insert(0, "__ROOT__/operators")
    from critical_eigvec_id import build_tau_np
    tau, ns = build_tau_np(64, 14)   # (N=64, c=14) -> real symmetric (2N+1)x(2N+1) matrix; ns = index grid
tau IS A+B-P on the test basis; RH <=> tau stays PSD for ALL (N,c). It is Arb-certified PD at finite (N,c).
MEASURED CRUX (2026-07-19, results/margin_scaling.json): the certified positivity margin (min LDL^T pivot)
PLATEAUS as N grows at fixed c (tame — a uniform lower bound in the test-basis limit), but COLLAPSES as the
prime cutoff c grows (razor-thin: at N=64 the margin falls from 4e-13 at c=8 to 8e-22 at c=14 and needs
>1000-bit precision beyond c=25). So the REAL wall is the c->infinity PRIME-SIDE limit, not the N-limit.
Your job is a STRUCTURAL attack on THAT: explain/bound how the prime side P eats the archimedean margin as
primes accumulate; find a cohomological/intersection signature that FORCES positivity uniformly in c; find a
factorization tau=M^*M with M built from zeta/prime data that stays valid as c grows; find why the near-null
pivot appears where it does and whether it can ever cross 0. Print a concrete diagnostic tied to build_tau_np
across SEVERAL c values. Toy constructions disconnected from build_tau_np, or that only vary N, are rejected."""
PROMPT=PROMPT.replace("__ELIM__",elim).replace("__TRIED__",tried).replace("__ROOT__",ROOT)

def call_claude(p):
    env=dict(os.environ); env.pop("ANTHROPIC_API_KEY",None)
    env["CLAUDE_CONFIG_DIR"]=os.path.expanduser("~/.claude-cobo")
    for args in ([CLAUDE,"--exclude-dynamic-system-prompt-sections","-p",p],[CLAUDE,"-p",p]):
        try:
            r=subprocess.run(args,capture_output=True,text=True,timeout=800,env=env)
            if r.returncode==0 and r.stdout.strip(): return r.stdout
        except Exception: continue
    return None

def strip_ansi(t):
    """ollama's TUI injects ANSI cursor codes; strip them so the JSON parses."""
    t=re.sub(r'\x1b\[[0-9;]*[A-Za-z]','',t)   # CSI sequences
    t=re.sub(r'\x1b\][^\x07]*\x07','',t)      # OSC sequences
    return t.replace('\r','')

def call_free(p):
    """FREE routing: Olympus GPU Ollama HTTP API (qwen2.5-coder:32b). Zero cash, zero Max quota.
    PRIMARY generator — Claude is only the rare capped escalation.
    Uses the HTTP API (/api/generate), NOT `ollama run`: the interactive TUI wrap-redraws lines and
    duplicates fragments, which broke json.loads (956 dead cycles). The API returns clean JSON."""
    for model in ["qwen2.5-coder:32b","qwen3:32b"]:
        payload=json.dumps({"model":model,"prompt":p,"stream":False,
                            "options":{"num_predict":4096,"temperature":0.7}})
        try:
            if ON_MAC:  # pipe payload over ssh stdin to Olympus's local ollama (no shell-escaping)
                r=subprocess.run(["ssh","olympus","curl -s --max-time 600 localhost:11434/api/generate -d @-"],
                                 input=payload,capture_output=True,text=True,timeout=660)
                raw=r.stdout
            else:       # daemon home: hit localhost directly
                req=urllib.request.Request("http://localhost:11434/api/generate",
                        data=payload.encode(),headers={"Content-Type":"application/json"})
                raw=urllib.request.urlopen(req,timeout=660).read().decode()
            resp=json.loads(raw).get("response","")
            resp=re.sub(r'<think>.*?</think>','',resp,flags=re.S)  # qwen3 emits reasoning blocks
            if re.search(r'\[.*\]',resp,re.S): return resp
        except Exception: continue
    return None

# FREE first
out=call_free(PROMPT); paid=False
if not out or not re.search(r'\[.*\]',out,re.S):
    ok,why=budget.can_spend_paid()
    if not ok:
        print(f"[{now()}] free route empty AND paid blocked ({why}) — skip generation, loop stays FREE"); 
        raise SystemExit(0)
    out=call_claude(PROMPT); paid=True
    budget.record_paid(success=bool(out and re.search(r'\[.*\]',out,re.S)))
if not out:
    open(f"{ENG}/GENERATE_PROMPT_pending.txt","w").write(PROMPT)
    print(f"[{now()}] Claude CLI unavailable — prompt saved for next session"); raise SystemExit(0)
m=re.search(r'\[.*\]',out,re.S)
try: methods=json.loads(m.group(0) if m else out, strict=False)
except Exception as e:
    open(f"{ENG}/last_gen_raw.txt","w").write(out); print(f"[{now()}] parse fail ({e}); raw saved"); raise SystemExit(0)
def smoke(path):
    """Fast run: (ok, tail). Skips clogging the queue with code that crashes on sight."""
    try:
        r=subprocess.run(["python3",path],capture_output=True,text=True,timeout=60,cwd=ROOT)
        return r.returncode==0, (r.stderr or r.stdout)[-600:]
    except subprocess.TimeoutExpired:
        return True, "(runs >60s — kept; real cycle allows 20min)"  # slow != broken
    except Exception as e:
        return False, str(e)

def repair(code,tail):
    """One FREE repair pass — hand the crash back to the model. Zero quota."""
    fix=call_free(f"This python3 script crashed. Return ONLY the corrected COMPLETE script "
                  f"(no prose, no markdown fences), same intent, must run clean and PRINT a clear "
                  f"numeric/structural signal.\n\n=== SCRIPT ===\n{code}\n\n=== ERROR ===\n{tail}")
    if not fix: return None
    fix=re.sub(r'^```[a-z]*\n?|```$','',fix.strip(),flags=re.M)
    return fix if "def " in fix or "import " in fix else None

added=0
for M in methods:
    nm=re.sub(r'[^a-zA-Z0-9_]','',M.get("name","auto_x"))[:40]
    code=M.get("python_code","")
    if not code or any(x["name"]==nm for x in reg["methods"]): continue
    if "build_tau_np" not in code:  # must attack the REAL operator, not a toy matrix
        print(f"[{now()}] x {nm} skipped (does not touch build_tau_np — toy)"); continue
    p=f"{ROOT}/operators/{nm}.py"; open(p,"w").write(code)
    ok,tail=smoke(p)
    if not ok:  # one free repair attempt before giving up on it
        fixed=repair(code,tail)
        if fixed:
            open(p,"w").write(fixed); ok,tail=smoke(p)
    if not ok:
        print(f"[{now()}] x {nm} skipped (crashes, unrepaired): {tail[:80]}"); continue
    reg["methods"].append({"name":nm,"code":f"operators/{nm}.py","venv":"sys","family":M.get("family","auto"),
        "status":"pending","yield":0.8,"born":now(),"rationale":M.get("rationale","")[:200],"novelty":M.get("novelty","")[:200]})
    added+=1
    print(f"[{now()}] + {nm} ({M.get('family')}): {M.get('rationale','')[:100]}")
json.dump(reg,open(REG,"w"),indent=2)
print(f"[{now()}] generation via {'CLAUDE(paid)' if paid else 'FREE'} added {added}; live queue now {sum(1 for x in reg['methods'] if x['status'] in ('pending','lead'))}")
