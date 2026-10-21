#!/usr/bin/env python3
"""Несколько решений × наборы тестов через настоящий интерактор (Windows), таблица средних KPI по семействам.

  python scripts/multi_bench.py exe1 exe2 ... -- dir1 dir2 ...  [-j 8]
Каталог может содержать *.in или подкаталоги-семейства с *.in. Результаты по тестам — build/multi_last.json.
TL: CPU > 5 с отмечается в столбце maxcpu (KPI не обнуляется — смотреть отдельно)."""
import json, os, sys
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from win_bench import run_one

args = sys.argv[1:]
j = 8
if '-j' in args:
    k = args.index('-j'); j = int(args[k + 1]); del args[k:k + 2]
i = args.index('--')
exes = [os.path.abspath(e) for e in args[:i]]
files = []
for d in args[i + 1:]:
    d = os.path.abspath(d)
    direct = sorted(f for f in os.listdir(d) if f.endswith('.in'))
    for f in direct: files.append((os.path.basename(d), os.path.join(d, f)))
    for sub in sorted(os.listdir(d)):
        p = os.path.join(d, sub)
        if os.path.isdir(p):
            for f in sorted(os.listdir(p)):
                if f.endswith('.in'): files.append((sub, os.path.join(p, f)))
jobs = [(e, f) for e in exes for _, f in files]
with ThreadPoolExecutor(j) as ex:
    res = dict(zip(jobs, ex.map(lambda ef: run_one(ef[0], ef[1], 60), jobs)))
json.dump({f"{os.path.basename(e)}|{f}": res[(e, f)] for e, f in jobs}, open('build/multi_last.json', 'w'))
fams = []
for fam, _ in files:
    if fam not in fams: fams.append(fam)
names = [os.path.basename(e)[:-4][:10] for e in exes]
print('family'.ljust(14), 'n'.rjust(3), ''.join(f'{n:>11}' for n in names))
tot = [0.0] * len(exes); cnt = 0
for fam in fams:
    fs = [f for fa, f in files if fa == fam]
    m = [sum(res[(e, f)][0] for f in fs) / len(fs) for e in exes]
    for k, x in enumerate(m): tot[k] += x * len(fs)
    cnt += len(fs)
    print(fam.ljust(14), str(len(fs)).rjust(3), ''.join(f'{x:11.4f}' for x in m))
print('ALL'.ljust(14), str(cnt).rjust(3), ''.join(f'{x / cnt:11.4f}' for x in tot))
print('maxcpu'.ljust(18), ''.join(f'{max(res[(e, f)][1] for _, f in files):11.2f}' for e in exes))
print('fails'.ljust(18), ''.join(f'{sum(1 for _, f in files if res[(e, f)][0] == 0):11d}' for e in exes))
