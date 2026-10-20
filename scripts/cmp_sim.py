#!/usr/bin/env python3
"""Потестовое сравнение двух выводов стенда sim2.exe.

  build\\sim_0.exe materials\\tests\\*.in > build\\a.txt
  build\\sim_1.exe materials\\tests\\*.in > build\\b.txt
  python scripts/cmp_sim.py build\\a.txt build\\b.txt        # A = моё, B = эталон
"""
import os
import sys


def load(path):
    d = {}
    for line in open(path):
        f = line.split()
        if len(f) != 7:
            continue
        name = os.path.basename(f[0])
        key = os.path.splitext(name)[0]
        d[key] = dict(w=float(f[1]), kpi=float(f[2]), T=int(f[3]), digs=int(f[4]), ok=int(f[6]))
    return d


def sort_key(k):
    return (0, int(k)) if k.isdigit() else (1, k)


def main():
    a = load(sys.argv[1])
    b = load(sys.argv[2])
    keys = sorted(set(a) & set(b), key=sort_key)
    print(f"{'test':>22} {'w':>4} {'A':>8} {'B':>8} {'d*w':>8} {'TA':>8} {'TB':>8} {'dT':>6}")
    tot = 0.0
    for k in keys:
        ra, rb = a[k], b[k]
        d = (ra["kpi"] - rb["kpi"]) * ra["w"]
        tot += d
        flag = "" if ra["ok"] else " WA!"
        print(f"{k:>22} {ra['w']:>4g} {ra['kpi']:>8.5f} {rb['kpi']:>8.5f} {d:>8.4f} "
              f"{ra['T']:>8} {rb['T']:>8} {ra['T'] - rb['T']:>6}{flag}")
    print(f"\nвзвешенная разница A−B: {tot:+.4f}")


if __name__ == "__main__":
    main()
