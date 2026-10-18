#!/usr/bin/env python3
"""Портфель против bf-L на build/fuzz2 (360 мелких случайных тестов): средняя разница и худшие случаи. Доп. аргументы — для sim (например --pset ...)."""
import subprocess, os, re, glob, sys
exe = os.path.abspath('build/bin/sim.exe')
extra = sys.argv[1:]
rows = []
for sd in sorted(glob.glob('build/fuzz2/*')):
    names = sorted([os.path.basename(f) for f in glob.glob(sd + '/*.in')], key=lambda x: (len(os.path.splitext(x)[0]), os.path.splitext(x)[0]))
    def run(args):
        o = subprocess.run([exe] + args + ['-v', sd], capture_output=True, text=True).stdout
        return [float(x) for x in re.findall(r'kpi=([0-9.]+)', o)]
    bf = run(['--pol', 'bf-L']); pt = run(['--port'] + extra)
    for n, a, b in zip(names, bf, pt): rows.append((b - a, sd, n, a, b))
rows.sort()
print("mean diff %.4f over %d tests; worst:" % (sum(r[0] for r in rows) / len(rows), len(rows)))
for r in rows[:5]: print("  %+.4f %s %s bf=%.4f port=%.4f" % r)
