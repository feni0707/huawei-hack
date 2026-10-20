#!/usr/bin/env python3
"""runroll.py "<ENV k=v k=v>;<ENV2>" file...  : прогон build/bin/lab_roll.exe с разными окружениями, печать T/kpi по файлам"""
import sys, os, subprocess, concurrent.futures as cf
exe = os.path.abspath('build/bin/lab_roll.exe')
cfgs = [c.strip() for c in sys.argv[1].split(';')]
files = sys.argv[2:]
def run(cfg, f):
    env = dict(os.environ)
    for k in ('ROLLH', 'ROLLMODE', 'ROLLK', 'ROLLMIN', 'ACCMIN', 'PSHARE', 'PCNT', 'ROLLM', 'HACCMIN'):
        env.pop(k, None)
    for kv in cfg.split():
        k, v = kv.split('=', 1); env[k] = v
    out = subprocess.run([exe, f], capture_output=True, text=True, env=env).stdout
    for l in out.splitlines():
        if 'kpi=' in l:
            t = l.split(); return t[2][4:], t[3][2:]
    return '?', '?'
with cf.ThreadPoolExecutor(8) as ex:
    futs = {(c, f): ex.submit(run, c, f) for c in cfgs for f in files}
for f in files:
    print(os.path.basename(f).ljust(10), ' | '.join(f"{(c or 'base')[:34]:34} kpi={futs[(c, f)].result()[0]} T={futs[(c, f)].result()[1]}" for c in cfgs))

