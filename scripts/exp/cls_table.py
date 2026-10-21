#!/usr/bin/env python3
"""cls_table.py t thr bucket cap out [--pre11] : таблица классов "pid класс" для теста t (как в scorefit_ov_build.py):
класс = 1 + смерть // bucket (не больше cap), если труба живёт >= thr команд (bucket 0 — просто 1); выжившие — 0.
--pre11: трубы общего с 11 префикса (первые 4999 команд) — по правилу 11 (50, 0), как в посылке (биты самого короткого теста)."""
import sys
def load(t):
    tok = open(f'materials/tests/{t}.in').read().split(); n = int(tok[1]); return list(map(int, tok[2:2 + n]))
def cls(c, thr, bucket, cap):
    birth = []; death = {}
    for q, x in enumerate(c):
        if x > 0: birth.append(q)
        else: death[-x] = q
    out = {}
    for pid in range(1, len(birth) + 1):
        d = death.get(pid)
        if d is not None and d - birth[pid - 1] >= thr: out[pid] = min(1 + d // bucket, cap) if bucket else 1
    return out
t, thr, bucket, cap, outp = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
c = load(t); k = cls(c, thr, bucket, cap)
if '--pre11' in sys.argv:
    c11 = load(11); pre = sum(1 for x in c[:4999] if x > 0) + 1
    k11 = cls(c11, 50, 0, 1)
    k = {p: v for p, v in k.items() if p >= pre}; k.update({p: v for p, v in k11.items() if p < pre})
with open(outp, 'w') as f:
    for p in sorted(k): f.write(f'{p} {k[p]}\n')
