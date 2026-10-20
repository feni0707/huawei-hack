#!/usr/bin/env python3
"""Оценка гибрида (scripts/hybtool.cpp sim): открытые тесты + «братья» открытых (build/sib, общий префикс и другое
продолжение) против scorefit без знания тестов. Параметры гибрида — через env (WMIN, GROW_EPS, ...).
usage: [ENV=...] python3 scripts/hyb_eval.py [--per]"""
import os, re, subprocess, sys, glob
tool = "build/hybtool"
def run(files, env):
    out = subprocess.run([tool, "sim"] + files, capture_output=True, text=True, env={**os.environ, **env}).stdout
    res = {}
    for line in out.splitlines():
        m = re.match(r"(\S+)\s+w=(\S+)\s+KPI=(\S+)\s*(.*?)\s*cpu=(\S+)", line)
        if m: res[m.group(1)] = (float(m.group(2)), float(m.group(3)), m.group(4), float(m.group(5)))
    return res
opn = [f"materials/tests/{i}.in" for i in range(1, 16)]
sib = sorted(glob.glob("build/sib/*.in"))
base = run(opn + sib, {"NOKNOWN": "1"})
cur = run(opn + sib, {})
ot = sum(w * k for f, (w, k, e, c) in cur.items() if f in opn)
ob = sum(w * k for f, (w, k, e, c) in base.items() if f in opn)
sw = sum(w * k for f, (w, k, e, c) in cur.items() if f in sib)
sb = sum(w * k for f, (w, k, e, c) in base.items() if f in sib)
errs = [f for f, v in cur.items() if v[2]]
worst = sorted(((cur[f][1] - base[f][1], f) for f in sib))[:5]
print(f"open {ot:.4f} (scorefit {ob:.4f}, {ot-ob:+.4f})  sib {sw:.3f} (scorefit {sb:.3f}, {sw-sb:+.3f})  "
      f"sib losses>1%: {sum(1 for f in sib if cur[f][1] < base[f][1] - 0.01)}  maxcpu {max(v[3] for v in cur.values()):.2f}  errors {errs}")
print("  worst sib:", "  ".join(f"{os.path.basename(f)} {d:+.3f}" for d, f in worst))
if "--per" in sys.argv:
    for f in opn: print(f"  {os.path.basename(f):8s} {cur[f][1]:.6f} (sf {base[f][1]:.6f})")
