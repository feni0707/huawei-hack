#!/usr/bin/env python3
"""rollcmp.py "<cfgA>;<cfgB>;..." glob... : lab_roll на всех файлах (параллельно), сводка по подпапкам и по файлам, где различие
 (по умолчанию первая конфигурация — база). Печатает среднее kpi по группе (родительская папка) и список различий."""
import sys, os, subprocess, glob, collections, concurrent.futures as cf
exe = os.path.abspath('build/bin/lab_roll.exe')
cfgs = [c.strip() for c in sys.argv[1].split(';')]
files = sorted(sum([glob.glob(g, recursive=True) for g in sys.argv[2:]], []))
def run(cfg, f):
    env = dict(os.environ)
    for k in ('ROLLH', 'ROLLMODE', 'ROLLK', 'ROLLMIN', 'ACCMIN', 'PSHARE', 'PCNT', 'ROLLM', 'HACCMIN'):
        env.pop(k, None)
    for kv in cfg.split():
        k, v = kv.split('=', 1); env[k] = v
    out = subprocess.run([exe, f], capture_output=True, text=True, env=env).stdout
    for l in out.splitlines():
        if 'kpi=' in l:
            return float(l.split()[2][4:])
    return 0.0
with cf.ThreadPoolExecutor(12) as ex:
    futs = {(c, f): ex.submit(run, c, f) for c in cfgs for f in files}
R = {k: v.result() for k, v in futs.items()}
grp = collections.defaultdict(list)
for f in files: grp[os.path.basename(os.path.dirname(f))].append(f)
print('group'.ljust(16), 'n'.rjust(3), ''.join(f'{(c or "base")[:22]:>24}' for c in cfgs))
tot = [0.0] * len(cfgs)
for g, fs in sorted(grp.items()):
    m = [sum(R[(c, f)] for f in fs) / len(fs) for c in cfgs]
    print(g.ljust(16), str(len(fs)).rjust(3), ''.join(f'{x:24.5f}' for x in m))
    for i, x in enumerate(m): tot[i] += x * len(fs)
print('ALL'.ljust(16), str(len(files)).rjust(3), ''.join(f'{x / len(files):24.5f}' for x in tot))
for f in files:
    b = R[(cfgs[0], f)]
    d = [R[(c, f)] - b for c in cfgs[1:]]
    if any(abs(x) > 1e-5 for x in d):
        print('  diff', os.path.relpath(f).ljust(34), f'base={b:.5f}', ' '.join(f'{x:+.5f}' for x in d))

