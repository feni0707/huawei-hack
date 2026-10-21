#!/usr/bin/env python3
"""Парное сравнение двух конфигураций варианта A (build/lab_a) на большом наборе тестов.

Для каждого класса: средний KPI обеих конфигураций, средняя разность (B − A), стандартная ошибка,
число тестов, где B лучше / хуже / равно. Нужен, чтобы отличить реальное улучшение от шума:
изменение значимо, если |разность| > 2·SE и wins/losses сильно несимметричны.

  python3 scripts/paired_a.py "pol=0" "pol=0 tieL=1 thr=7"
  python3 scripts/paired_a.py --dirs "build/gen_s*" "build/gen_big/s*" "pol=0" "pol=6 m=5 p=2"
"""
import argparse
import glob
import math
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LAB = os.path.join(ROOT, "build", "lab_a")
CLASSES = ["alternating", "periodic", "phases", "random_bursts", "sorted_groups", "short_long"]


def run(params, files):
    out = subprocess.run([LAB] + params + files, capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit(out.stderr)
    return {l.split()[0]: float(l.split()[1]) for l in out.stdout.strip().split("\n") if l}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--dirs", nargs="*", default=["build/gen_s*", "build/gen_big/s*"])
    ap.add_argument("--open", action="store_true", help="добавить открытые тесты отдельной строкой")
    args = ap.parse_args()
    pa, pb = args.a.split(), args.b.split()
    print(f"A = {args.a!r}\nB = {args.b!r}")
    print(f"{'class':14s} {'n':>4} {'A':>8} {'B':>8} {'B-A':>9} {'SE':>8} {'B>A':>4} {'B<A':>4} {'B=A':>4}")
    for c in CLASSES:
        files = []
        for d in args.dirs:
            files += sorted(glob.glob(os.path.join(ROOT, d, c, "*.in")))
        if not files:
            continue
        ka, kb = run(pa, files), run(pb, files)
        diffs = [kb[f] - ka[f] for f in files]
        n = len(diffs)
        mean = sum(diffs) / n
        var = sum((d - mean) ** 2 for d in diffs) / max(1, n - 1)
        se = math.sqrt(var / n)
        w = sum(d > 1e-9 for d in diffs)
        l = sum(d < -1e-9 for d in diffs)
        print(f"{c:14s} {n:4d} {sum(ka.values())/n:8.4f} {sum(kb.values())/n:8.4f} {mean:+9.5f} {se:8.5f} {w:4d} {l:4d} {n-w-l:4d}")
    if args.open:
        files = sorted(glob.glob(os.path.join(ROOT, "materials", "tests", "*.in")), key=lambda p: int(os.path.basename(p)[:-3]))
        ka, kb = run(pa, files), run(pb, files)
        sa = sum(float(open(f).readline().split()[0]) * ka[f] for f in files)
        sb = sum(float(open(f).readline().split()[0]) * kb[f] for f in files)
        print(f"{'open Σw':14s} {len(files):4d} {sa:8.4f} {sb:8.4f} {sb - sa:+9.5f}")


if __name__ == "__main__":
    main()
