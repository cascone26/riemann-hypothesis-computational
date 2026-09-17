#!/usr/bin/env python3
"""build_rh_visualizer.py -- generates rh_ledger_visualizer.html, a single self-contained
canvas-based pan+zoom map of the RH ledger's 1579 sections + their cross-references.

Design lifted from two of Tony's recent builds:
  - SCALAR's canvas-over-SVG fix (2026-08-20): thousands of nodes/edges as batched canvas
    draws, not DOM/SVG elements -- SVG choked at ~67k elements, canvas did the same job
    at ~13ms. This ledger has 1579 nodes + ~6500 edges; canvas is the right call from the
    start, not a fix applied after it's already slow.
  - SCALAR's "no roll, just pan + continuous zoom" navigation model (renormalization):
    scroll to zoom (scale = the depth axis), drag to pan -- no 3D rotation, because the
    thing being explored (citation density across 1579 sections) has no "orientation",
    just more or less detail at different scales. Level-of-detail is scale-gated: dots
    only, at low zoom; titles + citation edges fade in above a zoom threshold.

Regenerate after re-running rh_ledger_graph.py build (this script reads its JSON output).

Usage: python3 tools/build_rh_visualizer.py
"""
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))
GRAPH_JSON = os.path.join(HERE, ".rh_ledger_graph.json")
OUT_HTML = os.path.join(HERE, "rh_ledger_visualizer.html")

with open(GRAPH_JSON, encoding="utf-8") as fh:
    graph = json.load(fh)

sections = graph["sections"]
# stable layout: order by section number, wrap into a grid so the map has real 2D
# extent to pan around in (a single 1579-wide line would make "pan" meaningless).
COLS = 40
nodes = []
id_to_index = {}
for i, s in enumerate(sections):
    id_to_index[s["id"]] = i
    col = i % COLS
    row = i // COLS
    status = "frontier" if (s["has_frontier_language"] and not s["has_closed_language"]) else (
        "lean" if s["lean_files"] else ("closed" if s["has_closed_language"] else "plain")
    )
    nodes.append({"id": s["id"], "col": col, "row": row, "title": s["title"][:140], "status": status, "line": s["line"]})

edges = []
for s in sections:
    a = id_to_index.get(s["id"])
    for c in s["cites"]:
        b = id_to_index.get(c)
        if a is not None and b is not None and a != b:
            edges.append([a, b])

DATA = {"nodes": nodes, "edges": edges, "cols": COLS}

HTML = """<!doctype html><html><head><meta charset="utf-8"><title>RH Ledger Map</title>
<style>
  html,body{margin:0;height:100%;background:#0b0d12;overflow:hidden;font-family:-apple-system,Helvetica,Arial,sans-serif;}
  canvas{display:block;cursor:grab;}
  #hud{position:fixed;top:10px;left:10px;color:#cfd6e4;font-size:12px;background:rgba(20,22,30,.7);
       padding:8px 10px;border-radius:8px;line-height:1.5;pointer-events:none;}
  #panel{position:fixed;top:10px;right:10px;max-width:360px;color:#eef1f7;font-size:13px;
         background:rgba(20,22,30,.9);padding:10px 12px;border-radius:8px;display:none;
         border:1px solid rgba(255,255,255,.08);}
  #panel b{color:#8ecbff;}
  .legend span{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:5px;vertical-align:middle;}
</style></head>
<body>
<canvas id="c"></canvas>
<div id="hud">
  RH Ledger Map -- __N_NODES__ sections, __N_EDGES__ citations<br>
  scroll = zoom (renormalize) &nbsp; drag = pan &nbsp; click a node = details<br>
  <span class="legend"><span style="background:#5ee6a8"></span>closed</span>
  <span class="legend"><span style="background:#7db8ff"></span>lean-backed</span>
  <span class="legend"><span style="background:#ffb15e"></span>frontier</span>
  <span class="legend"><span style="background:#5a6070"></span>plain</span>
</div>
<div id="panel"></div>
<script>
const DATA = __DATA_JSON__;
const canvas = document.getElementById('c');
const ctx = canvas.getContext('2d');
const panel = document.getElementById('panel');
let W, H;
function fitCanvas(){ W = canvas.width = innerWidth; H = canvas.height = innerHeight; }
fitCanvas();  // dimensions only, no draw yet -- draw() isn't defined and scale/offX
              // aren't initialized this early; calling it here is a TDZ crash
              // (the exact bug class BLE's own postmortem flagged -- see
              // [[filament-viewer-project]] visible-error-banner note).

const CELL = 26;
let scale = 1, offX = W/2 - (DATA.cols*CELL)/2, offY = 60;
let dragging = false, lastX=0, lastY=0;

const COLORS = {closed:'#5ee6a8', lean:'#7db8ff', frontier:'#ffb15e', plain:'#5a6070'};

function worldToScreen(x,y){ return [x*scale + offX, y*scale + offY]; }

function draw(){
  ctx.fillStyle = '#0b0d12';
  ctx.fillRect(0,0,W,H);

  const showLabels = scale > 5;
  const showEdges = scale > 2.2;

  if (showEdges) {
    ctx.strokeStyle = 'rgba(140,160,200,0.10)';
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (const [a,b] of DATA.edges){
      const na = DATA.nodes[a], nb = DATA.nodes[b];
      const [ax,ay] = worldToScreen(na.col*CELL+CELL/2, na.row*CELL+CELL/2);
      const [bx,by] = worldToScreen(nb.col*CELL+CELL/2, nb.row*CELL+CELL/2);
      if (ax<-50||ax>W+50||bx<-50||bx>W+50) continue;
      ctx.moveTo(ax,ay); ctx.lineTo(bx,by);
    }
    ctx.stroke();
  }

  const r = Math.max(2, Math.min(9, scale*3));
  for (const n of DATA.nodes){
    const [x,y] = worldToScreen(n.col*CELL+CELL/2, n.row*CELL+CELL/2);
    if (x<-20||x>W+20||y<-20||y>H+20) continue;
    ctx.fillStyle = COLORS[n.status] || COLORS.plain;
    ctx.beginPath(); ctx.arc(x,y,r,0,Math.PI*2); ctx.fill();
    if (showLabels){
      ctx.fillStyle = 'rgba(230,235,245,0.85)';
      ctx.font = '10px -apple-system,Helvetica,Arial';
      ctx.fillText(n.id, x+r+3, y+3);
    }
  }
}

canvas.addEventListener('wheel', e => {
  e.preventDefault();
  const factor = Math.exp(-e.deltaY * 0.001);
  const mx = e.clientX, my = e.clientY;
  const wx = (mx - offX) / scale, wy = (my - offY) / scale;
  scale = Math.max(0.15, Math.min(40, scale * factor));
  offX = mx - wx*scale; offY = my - wy*scale;
  draw();
}, {passive:false});

canvas.addEventListener('mousedown', e => { dragging=true; lastX=e.clientX; lastY=e.clientY; canvas.style.cursor='grabbing'; });
addEventListener('mouseup', () => { dragging=false; canvas.style.cursor='grab'; });
addEventListener('mousemove', e => {
  if(!dragging) return;
  offX += e.clientX-lastX; offY += e.clientY-lastY;
  lastX=e.clientX; lastY=e.clientY;
  draw();
});
canvas.addEventListener('click', e => {
  for (const n of DATA.nodes){
    const [x,y] = worldToScreen(n.col*CELL+CELL/2, n.row*CELL+CELL/2);
    if (Math.hypot(e.clientX-x, e.clientY-y) < 10){
      panel.style.display='block';
      panel.innerHTML = '<b>'+n.id+'</b> (line '+n.line+', status: '+n.status+')<br>'+n.title;
      return;
    }
  }
  panel.style.display='none';
});

function onResize(){ fitCanvas(); draw(); }
addEventListener('resize', onResize);
draw();  // initial paint, now that scale/offX/offY/draw all exist
</script>
</body></html>"""

html = HTML.replace("__DATA_JSON__", json.dumps(DATA)).replace("__N_NODES__", str(len(nodes))).replace("__N_EDGES__", str(len(edges)))
with open(OUT_HTML, "w", encoding="utf-8") as fh:
    fh.write(html)
print(f"wrote {OUT_HTML} ({len(nodes)} nodes, {len(edges)} edges, {len(html)} bytes)")
