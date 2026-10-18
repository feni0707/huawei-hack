#!/usr/bin/env python3
"""Случайный поиск параметров SCORED-политики (sc) на наборах классов. Печатает лучшие.
  python scripts/rsearch.py [N]"""
import os, random, subprocess, sys, re
N = int(sys.argv[1]) if len(sys.argv) > 1 else 100
SETS = ["materials/tests"]
for s in (1, 2, 3):
    SETS += [f"build/gen2_s{s}/twophase", f"build/gen2_s{s}/invert", f"build/gen_s{s}/alternating", f"build/gen_s{s}/phases"]
for s in (1, 2):
    SETS += [f"build/gen3_s{s}/sawtooth", f"build/gen3_s{s}/fifo", f"build/gen3_s{s}/desc"]
W = {"tests": 0.3 / 32.1, "twophase": 0.3, "invert": 0.1, "alternating": 0.15, "phases": 0.05, "sawtooth": 0.05, "fifo": 0.05, "desc": 0.05}
def run(pset):
    cmd = [os.path.abspath("build/bin/sim_rs.exe"), "--pol", "sc", "--pool", "build/pool_sc.txt", "--line"] + (["--pset", "sc:" + pset] if pset else []) + SETS
    out = subprocess.run(cmd, capture_output=True, text=True).stdout
    d = {m.group(1): float(m.group(2)) for m in re.finditer(r"(\w+)=([0-9.]+)", out)}
    return sum(W[k] * d[k] for k in W), d
base_bf = run("")  # sc с нулевыми весами = best fit-подобный
print("sc default:", round(base_bf[0], 5), base_bf[1])
res = []
rng = random.Random(7)
for i in range(N):
    th = rng.choice([2, 3, 4, 6, 8, 10, 12])
    sw0 = rng.choice([0, 20, 100, 400, 2000])
    swsmall = rng.choice([0, 5, 20, 100, 400])
    swsize = rng.choice([-30, -10, -3, -1, 0, 1, 3, 10, 30])
    ps = f"th={th},sw0={sw0},swsmall={swsmall},swsize={swsize}"
    sc, d = run(ps)
    res.append((sc, ps, d))
    if i % 20 == 19: print(i + 1, "best so far", max(res)[0], max(res)[1], flush=True)
res.sort(reverse=True)
for sc, ps, d in res[:8]:
    print(round(sc, 5), ps, {k: round(v, 4) for k, v in d.items()})
