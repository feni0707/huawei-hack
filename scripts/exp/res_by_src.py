#!/usr/bin/env python3
"""res_by_src.py base.exe new.exe [...] : средние KPI остаточных тестов по исходному открытому тесту (самый длинный общий
префикс) из build/multi_last.json (scripts/multi_bench.py). Печатает для каждого решения среднее и число хуже/лучше base."""
import json, os, sys
from collections import defaultdict
ex = sys.argv[1:]
d = json.load(open('build/multi_last.json'))
r = {e: {k.split('|', 1)[1]: x[0] for k, x in d.items() if k.startswith(e + '|')} for e in ex}
def load(fn):
    tok = open(fn).read().split(); n = int(tok[1]); return list(map(int, tok[2:2 + n]))
op = {t: load(f'materials/tests/{t}.in') for t in range(1, 16)}
g = defaultdict(list)
for f in r[ex[0]]:
    if not all(f in r[e] for e in ex): continue
    c = load(f); best = (0, 0)
    for t, o in op.items():
        k = 0; m = min(len(o), len(c))
        while k < m and o[k] == c[k]: k += 1
        if k > best[0]: best = (k, t)
    g[best[1]].append(tuple(r[e][f] for e in ex))
for t in sorted(g):
    v = g[t]; n = len(v)
    s = f'src {t:2d}: n {n:3d}  {ex[0]} {sum(x[0] for x in v) / n:.5f}'
    for i in range(1, len(ex)):
        s += f' | {ex[i]} {sum(x[i] for x in v) / n:.5f} ({sum(x[i] - x[0] for x in v) / n:+.5f}; -{sum(1 for x in v if x[i] < x[0] - 1e-9)} +{sum(1 for x in v if x[i] > x[0] + 1e-9)})'
    print(s)
