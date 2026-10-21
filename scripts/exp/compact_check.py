#!/usr/bin/env python3
"""Насколько вшитая раскладка (build/off/<t>.res, строка 4 — D по PID) раздувает траншею по ходу теста:
T_emb(t) против онлайн-границы b(t) = 20·ceil(max живых до t / 20). Печатает максимум и среднее превышение (доля b)."""
import sys
def load(t):
    f = open(f'materials/tests/{t}.in').read().split(); n = int(f[1]); return list(map(int, f[2:2 + n]))
for t in map(int, sys.argv[1:] or range(5, 16)):
    c = load(t)
    lines = open(sys.argv[0].replace('compact_check.py', '') + f'../../build/off/{t}.res').read().split('\n') if False else open(__import__('os').environ.get('RESDIR','build/off')+f'/{t}.res').read().split('\n')
    D = [0] + list(map(int, lines[3].split()))
    ln = [0]; live = 0; pk = 0; mx = 0; worst = 0; wt = 0; s = 0; cnt = 0; pid = 0
    for i, x in enumerate(c):
        if x > 0:
            pid += 1; ln.append(x); live += x; pk = max(pk, live); mx = max(mx, D[pid] + x)
        else: live -= ln[-x]
        b = (pk + 19) // 20 * 20; T = (mx + 19) // 20 * 20
        ex = (T - b) / max(b, 1)
        if ex > worst: worst = ex; wt = i
        s += ex; cnt += 1
    print(f'test {t:2d}: max excess {100*worst:7.1f}% of bound (at cmd {wt}/{len(c)}), mean {100*s/cnt:6.1f}%')
