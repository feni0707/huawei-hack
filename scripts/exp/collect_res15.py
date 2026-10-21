#!/usr/bin/env python3
"""Собрать остаточные тесты корпуса, начинающиеся с открытого теста T (общий префикс >= MINP), в build/resT/r/."""
import glob, os, shutil, sys
T = int(sys.argv[1]) if len(sys.argv) > 1 else 15
MINP = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
o = open(f'materials/tests/{T}.in').read().split(); o = list(map(int, o[2:2 + int(o[1])]))
out = f'build/res{T}/r'; os.makedirs(out, exist_ok=True)
cnt = 0
for f in sorted(glob.glob('build/corpus/res*/residual/*.in')):
    c = open(f).read().split(); n = int(c[1]); c = list(map(int, c[2:2 + n]))
    k = 0
    while k < min(len(o), len(c)) and o[k] == c[k]: k += 1
    if k >= MINP:
        cnt += 1
        name = os.path.basename(os.path.dirname(os.path.dirname(f))) + '_' + os.path.basename(f)
        shutil.copy(f, os.path.join(out, name))
        print(f, 'prefix', k, 'n', n)
print(cnt, 'tests ->', out)
