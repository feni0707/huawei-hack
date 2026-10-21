#!/usr/bin/env python3
"""cf_collect.py out.txt glob... : параллельный прогон build/bin/cfdump.exe (SAMPLE, MAXOPT из env) по тестам, строки R в out.txt"""
import sys, os, glob, subprocess, concurrent.futures as cf
out = sys.argv[1]; files = sorted(sum([glob.glob(g) for g in sys.argv[2:]], []))
exe = os.path.abspath('build/bin/cfdump.exe')
def one(f):
    return subprocess.run([exe, f], capture_output=True, text=True).stdout
with cf.ThreadPoolExecutor(int(os.environ.get('J', '16'))) as ex, open(out, 'w') as o:
    for i, r in enumerate(ex.map(one, files)):
        o.write(r)
print('files', len(files))
