#!/usr/bin/env python3
"""Разбор потерь решения: сколько «лишних» 20-метровых докопок сделано против теоретического минимума.

Минимум: T_min = 20*ceil(peak_live/20), где peak_live — максимум суммарной длины одновременно живых труб.
Он не зависит от стратегии (КПИ считается только по финальному состоянию, траншея не укорачивается),
поэтому единственное, что оптимизирует решение, — число докопок. Строка отчёта показывает, сколько
докопок «лишние» и сколько баллов они стоят.

  python scripts/loss_report.py build/bench/ref_open.csv [--dir materials/tests]
"""
import argparse
import csv
import math
import os
import sys


def test_stats(path):
    t = open(path).read().split()
    w, n = float(t[0]), int(t[1])
    cmds = list(map(int, t[2:2 + n]))
    alive = {}
    pid = 0
    cur = peak = 0
    dels = 0
    lays = 0
    for c in cmds:
        if c > 0:
            pid += 1
            alive[pid] = c
            cur += c
            peak = max(peak, cur)
            lays += 1
        else:
            cur -= alive.pop(-c)
            dels += 1
    return dict(w=w, n=n, live=cur, peak=peak, tmin=20 * math.ceil(peak / 20),
                lays=lays, dels=dels)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--dir", default="materials/tests")
    ap.add_argument("--set", default="open")
    args = ap.parse_args()

    stats = {}
    for f in sorted(os.listdir(args.dir)):
        if f.endswith(".in"):
            stats[f] = test_stats(os.path.join(args.dir, f))

    rows = [r for r in csv.DictReader(open(args.csv)) if r["set"] == args.set]
    by_sol = {}
    for r in rows:
        by_sol.setdefault(r["solution"], []).append(r)

    print(f"{'test':>5} {'w':>4} {'n':>7} {'lays':>7} {'dels':>7} {'peak':>9} {'Tmin':>9} "
          f"{'kpi_ub':>7}", end="")
    sols = list(by_sol)
    for s in sols:
        print(f" | {s[:12]:>12} {'T':>8} {'extra':>6} {'blk':>4}", end="")
    print()

    total = {s: 0.0 for s in sols}
    for f in sorted(stats, key=lambda x: int(os.path.splitext(x)[0])):
        st = stats[f]
        print(f"{os.path.splitext(f)[0]:>5} {st['w']:>4g} {st['n']:>7} {st['lays']:>7} {st['dels']:>7} "
              f"{st['peak']:>9} {st['tmin']:>9} {st['live'] / st['tmin']:>7.4f}", end="")
        for s in sols:
            r = next((x for x in by_sol[s] if x["test"] == f), None)
            if not r:
                print(f" | {'-':>12} {'-':>8} {'-':>6} {'-':>4}", end="")
                continue
            T = int(r["trench"])
            kpi = float(r["kpi"])
            extra = T - st["tmin"]
            total[s] += (st["live"] / st["tmin"] - kpi) * st["w"]
            print(f" | {s[:12]:>12} {T:>8} {extra:>6} {extra // 20:>4}", end="")
        print()
    print("\nразрыв до UB по взвешенной сумме:")
    for s in sols:
        print(f"  {s:>14}: {total[s]:.4f}")


if __name__ == "__main__":
    main()
