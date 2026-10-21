#!/usr/bin/env python3
"""delstats.py dir_or_files... : статистика процесса удалений по тестам — доля удалений среди 8 последних той же длины (lifo8),
средний квантиль возраста удаляемой среди живых (q: 0 — самая свежая, 1 — самая старая), доля удалений «только что уложенной» (pop)."""
import sys, os, glob, bisect, collections
def files_of(args):
    out = []
    for a in args:
        if a.endswith('.in'): out.append(a); continue
        out += sorted(glob.glob(os.path.join(a, '*.in'))) + sorted(glob.glob(os.path.join(a, '*', '*.in')))
    return out
def stats(path):
    f = open(path).read().split(); n = int(f[1]); c = list(map(int, f[2:2 + n]))
    ln = [0]; bylen = collections.defaultdict(list); alive = []  # alive sorted pids
    pid = 0; D = 0; hit = 0; qs = 0.0; pops = 0; prevp = 0
    import bisect
    for x in c:
        if x > 0:
            pid += 1; ln.append(x); bylen[x].append(pid); alive.append(pid); prevp = pid
        else:
            p = -x; L = ln[p]; st = bylen[L]; j = bisect.bisect_left(st, p)
            if len(st) - 1 - j < 8: hit += 1
            st.pop(j)
            k = bisect.bisect_left(alive, p); na = len(alive)
            qs += (na - 1 - k) / max(1, na - 1)   # доля живых моложе удаляемой
            alive.pop(k)
            pops += (p == prevp); prevp = 0
            D += 1
    return D, hit / max(D, 1), qs / max(D, 1), pops / max(D, 1)
rows = collections.defaultdict(list)
for p in files_of(sys.argv[1:]):
    fam = os.path.basename(os.path.dirname(p))
    rows[fam].append((os.path.basename(p),) + stats(p))
for fam, rs in rows.items():
    if len(rs) > 3:
        m = lambda k: sum(r[k] for r in rs) / len(rs)
        print(f'{fam:16s} n={len(rs):3d} lifo8={m(2):.2f} q={m(3):.2f} pop={m(4):.2f}')
    else:
        for r in rs: print(f'{fam:16s} {r[0]:10s} D={r[1]:6d} lifo8={r[2]:.2f} q={r[3]:.2f} pop={r[4]:.2f}')
