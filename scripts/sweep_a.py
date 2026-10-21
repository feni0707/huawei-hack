#!/usr/bin/env python3
"""Перебор параметров варианта A на оффлайн-симуляторе (build/lab_a, см. solutions/variant_a.cpp).

Симулятор делает ровно то же, что решение + интерактор, но без пайпов: весь набор тестов за ~секунды.
Метрики на конфигурацию:
  open — Σ w·KPI на materials/tests (макс 32.1);
  alt/per/pha/rnd/srt/s_l — средний KPI по классам gen-наборов (build/gen_s*/<class>/*.in, все seed вместе);
  obj — взвешенная цель для выбора конфигурации: 10·s_l + 4·alt + 4·pha + 2·rnd + 1·per + 1·srt, делённая на 22
        (веса отражают оценку состава закрытых тестов: 27–30 ~ как 15 (класс short_long), 16–26 ~ как 12–14);
  worst — минимум obj по отдельным seed (устойчивость).

Использование:
  python3 scripts/sweep_a.py "pol=0,1,2 side=0,1"              # декартово произведение
  python3 scripts/sweep_a.py --base "pol=0 tie=0" "side=0,1"   # фиксированная часть + перебор
  python3 scripts/sweep_a.py --top 20 "polS=0,1,2 polL=0,1,2 thr=8,16"
"""
import argparse
import glob
import itertools
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LAB = os.path.join(ROOT, "build", "lab_a")
CLASSES = ["alternating", "periodic", "phases", "random_bursts", "sorted_groups", "short_long"]
SHORT = {"alternating": "alt", "periodic": "per", "phases": "pha", "random_bursts": "rnd", "sorted_groups": "srt", "short_long": "s_l"}
CW = {"alternating": 4, "periodic": 1, "phases": 4, "random_bursts": 2, "sorted_groups": 1, "short_long": 10}


def natural(p):
    b = os.path.basename(p)
    return (len(b), b)


def load_sets():
    open_files = sorted(glob.glob(os.path.join(ROOT, "materials", "tests", "*.in")),
                        key=lambda p: int(os.path.basename(p).split(".")[0]))
    weights = {}
    for p in open_files:
        with open(p) as f:
            weights[p] = float(f.readline().split()[0])
    gen = {}   # (class, seedname) -> files
    for d in sorted(glob.glob(os.path.join(ROOT, "build", "gen_s*"))):
        seed = os.path.basename(d)
        for c in CLASSES:
            fs = sorted(glob.glob(os.path.join(d, c, "*.in")), key=natural)
            if fs:
                gen[(c, seed)] = fs
    return open_files, weights, gen


OPEN, WEIGHTS, GEN = load_sets()
ALL = OPEN + [f for fs in GEN.values() for f in fs]


def run_config(params):
    out = subprocess.run([LAB] + params + ALL, capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(out.stderr)
    kpi = {}
    for line in out.stdout.split("\n"):
        if line:
            name, k = line.split()[:2]
            kpi[name] = float(k)
    res = {"open": sum(WEIGHTS[p] * kpi[p] for p in OPEN)}
    per_class = {}
    for c in CLASSES:
        vals = [kpi[f] for (cc, s), fs in GEN.items() if cc == c for f in fs]
        per_class[c] = sum(vals) / len(vals)
        res[SHORT[c]] = per_class[c]
    res["obj"] = sum(CW[c] * per_class[c] for c in CLASSES) / sum(CW.values())
    seeds = sorted({s for (_, s) in GEN})
    worst = 1e9
    for s in seeds:
        o = 0.0
        for c in CLASSES:
            fs = GEN.get((c, s), [])
            if fs:
                o += CW[c] * sum(kpi[f] for f in fs) / len(fs)
        worst = min(worst, o / sum(CW.values()))
    res["worst"] = worst
    return res


def parse_grid(spec):
    keys, vals = [], []
    for tok in spec.split():
        k, v = tok.split("=")
        keys.append(k)
        vals.append(v.split(","))
    for combo in itertools.product(*vals):
        yield [f"{k}={v}" for k, v in zip(keys, combo)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("grid", help='например "pol=0,1,2 side=0,1"')
    ap.add_argument("--base", default="", help="фиксированные параметры: 'k=v k=v'")
    ap.add_argument("--top", type=int, default=0, help="показать только лучшие N по obj")
    ap.add_argument("-j", type=int, default=2)
    args = ap.parse_args()
    base = args.base.split() if args.base else []
    configs = [base + g for g in parse_grid(args.grid)]
    print(f"{len(configs)} configs, {len(OPEN)} open + {len(ALL) - len(OPEN)} gen files", file=sys.stderr)
    with ThreadPoolExecutor(args.j) as ex:
        results = list(ex.map(run_config, configs))
    rows = sorted(zip(configs, results), key=lambda r: -r[1]["obj"])
    if args.top:
        rows = rows[:args.top]
    print(f"{'open':>8} {'alt':>6} {'per':>6} {'pha':>6} {'rnd':>6} {'srt':>6} {'s_l':>6} {'obj':>7} {'worst':>7}  config")
    for cfg, r in rows:
        print(f"{r['open']:8.4f} {r['alt']:6.4f} {r['per']:6.4f} {r['pha']:6.4f} {r['rnd']:6.4f} {r['srt']:6.4f} {r['s_l']:6.4f} "
              f"{r['obj']:7.4f} {r['worst']:7.4f}  {' '.join(cfg)}")


if __name__ == "__main__":
    main()
