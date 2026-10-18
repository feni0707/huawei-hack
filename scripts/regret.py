#!/usr/bin/env python3
"""Сравнение портфеля с чистыми политиками на тестах: regret.py <pool.csv> <port_verbose.txt> <pol1,pol2,...> [worst N]"""
import csv, re, sys, collections
pure = {}
order = collections.defaultdict(list)
for r in csv.DictReader(open(sys.argv[1])):
    pure[(r['set'], r['test'], r['pol'])] = float(r['kpi'])
    if r['pol'] == 'bf-L':
        order[(r['set'].split('/')[-1], r['test'])].append(r['set'])
pols = sys.argv[3].split(',')
res = []
for l in open(sys.argv[2]):
    m = re.match(r'\s+(\S+)\s+(\S+)\s+w=(\S+)\s+kpi=([0-9.]+)', l)
    if m: res.append((m.group(1), m.group(2), float(m.group(4))))
cnt = collections.Counter(); rows = []
for name, test, k in res:
    i = cnt[(name, test)]; cnt[(name, test)] += 1
    label = order[(name, test)][i]
    best = max(pure[(label, test, p)] for p in pols)
    bf = pure[(label, test, 'bf-L')]
    rows.append((k - best, k, bf, best, label, test))
rows.sort()
n = int(sys.argv[4]) if len(sys.argv) > 4 else 12
for d, k, bf, best, label, test in rows[:n]:
    print(f"{label:28s} {test:6s} port={k:.4f} bf={bf:.4f} oracle={best:.4f} regret={d:+.4f}")
bycls = collections.defaultdict(list)
for d, k, bf, best, label, test in rows: bycls[label.split('/')[-1]].append(d)
print("mean regret by class:", {c: round(sum(v)/len(v), 4) for c, v in bycls.items()})
