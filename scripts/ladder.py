#!/usr/bin/env python3
"""ladder.py : лестница от чистого best fit до scorefit в лаборатории (scorefit_lab -DSIM), включая веса чистой версии.
Печатает взвешенный open (Σ w·KPI по materials/tests) и obj по классам hold-out (cls-веса как у feni)."""
import os, sys, glob, subprocess, collections, concurrent.futures as cf
exe = os.path.abspath('build/bin/lab.exe')
CLEAN = '1,-38.8,-22.4,13.5,47,-203,107,-7.25,2.78,2.26,-8.96,23.4,1.05e+03,69,0,46.7,27.9,50,0,0,0,0,0,0,0,20,50'
BF = ','.join(['1'] + ['0'] * 26)
stages = [
    ('S0 best fit (левый край)', '10=0 80=0 8=0', BF),
    ('S1 + остаток к молодому соседу', '10=0 80=0 8=1', BF),
    ('S2 + потоки (чётность)', '10=2 80=0 8=1', BF),
    ('S3 + зоны по длинам', '10=2 80=1 8=1', BF),
    ('S4 + линейная оценка (веса scorefit)', '10=2 80=1 8=1', CLEAN),
    ('S5 = S4 с P8 по умолчанию', '', CLEAN),
]
open_files = sorted(glob.glob('materials/tests/*.in'), key=lambda f: int(os.path.basename(f)[:-3]))
hold = []
for d in sys.argv[1:]:
    for cd in sorted(glob.glob(os.path.join(d, '*'))):
        if os.path.isdir(cd):
            for f in sorted(glob.glob(os.path.join(cd, '*.in'))): hold.append((os.path.basename(cd), f))
W = {'two_phase': 10, 'two_phase_real': 10, 'alternating': 4, 'lifo_mix': 4, 'phases': 4, 'sorted_groups': 2, 'periodic': 2, 'random_bursts': 2}
OW = {1: .1, 12: 4, 13: 4, 14: 4, 15: 10}
def run(sf, sfw, files):
    env = dict(os.environ)
    env['SFW'] = sfw
    if sf: env['SF'] = sf
    else: env.pop('SF', None)
    chunks = [files[j::12] for j in range(12)]
    def one(ch):
        out = subprocess.run([exe] + ch, capture_output=True, text=True, env=env).stdout
        return {l.split()[0]: float(l.split()[2][4:]) for l in out.splitlines() if ' w=' in l and 'kpi=' in l}
    r = {}
    with cf.ThreadPoolExecutor(12) as ex:
        for x in ex.map(one, chunks): r.update(x)
    return r
for name, sf, sfw in stages:
    ro = run(sf, sfw, open_files)
    openscore = sum(OW.get(int(os.path.basename(f)[:-3]), 1) * ro[f] for f in open_files)
    line = f'{name:40s} open={openscore:7.3f}'
    if hold:
        rh = run(sf, sfw, [f for _, f in hold])
        cls = collections.defaultdict(list)
        for c, f in hold: cls[c].append(rh[f])
        obj = sum(W.get(c, 1) * sum(v) / len(v) for c, v in cls.items())
        line += f'  hold-obj={obj:7.3f}  ' + ' '.join(f'{c[:5]}={sum(v)/len(v):.3f}' for c, v in sorted(cls.items()))
    print(line, flush=True)
