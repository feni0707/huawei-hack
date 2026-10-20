#!/usr/bin/env python3
"""Быстрый стенд: собирает scripts/sim2.cpp с решениями и сравнивает метрики по наборам тестов.

Решение включается в стенд через #include, поэтому все решения должны определять struct Engine с методами
place(int L, ll& N, ll& D) / remove(int pid) — как solutions/universal.cpp.

  python scripts/simrun.py --sol solutions/universal.cpp                 # open + build/gen/*
  python scripts/simrun.py --sol solutions/universal.cpp --sol solutions/scorefit.cpp  # сравнение
  python scripts/simrun.py --sol ... --sets build/gen/two_phase_real build/gen/alternating
  python scripts/simrun.py --sol ... --K "3=1.5,7=0"                    # ручки эксперимента

Разбивка по наборам: набор = папка; печатается средний KPI, худший тест, сумма CPU и провалы (WA).
"""
import argparse
import os
import re
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def find_vcvars():
    for cand in (r"C:\Program Files\Microsoft Visual Studio\2022\Community",
                 r"C:\Program Files\Microsoft Visual Studio\2022\Professional",
                 r"C:\Program Files\Microsoft Visual Studio\2022\Enterprise"):
        p = os.path.join(cand, "VC", "Auxiliary", "Build", "vcvars64.bat")
        if os.path.exists(p):
            return p
    sys.exit("vcvars64.bat not found")


def build(sol, out, defines=""):
    src = os.path.abspath(os.path.join(ROOT, "scripts", "sim2.cpp"))
    solpath = os.path.abspath(sol)
    # решения вида scorefit.cpp определяют struct Alloc вместо struct Engine
    text = open(solpath, encoding="utf-8", errors="replace").read()
    if "struct Alloc" in text and "struct Engine" not in text:
        defines += " /DUSE_ALLOC"
    # через .bat: cmd.exe искажает экранированные кавычки в /D...="путь"
    bat = os.path.join(ROOT, "build", f"_sim_{abs(hash((sol, out))) % 100000}.bat")
    with open(bat, "w", encoding="ascii") as f:
        f.write(f'@echo off\ncall "{find_vcvars()}" >nul 2>nul\ncd /d "{ROOT}"\n'
                f'cl /nologo /O2 /EHsc /std:c++17 /DSIM /Ibuild\\ref\\include '
                f'"/DUNIVSRC=\\"{solpath}\\"" {defines} /Fo:build\\ /Fe:{out} "{src}"\n')
    r = subprocess.run(["cmd", "/c", bat], capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(os.path.join(ROOT, out)):
        print(r.stdout[-4000:], r.stderr[-4000:])
        sys.exit(f"build failed for {sol}")


def tests_in(d):
    files = [f for f in os.listdir(d) if f.endswith(".in")]
    files.sort(key=lambda f: (0, int(os.path.splitext(f)[0])) if os.path.splitext(f)[0].isdigit()
               else (1, 0, f))
    return [os.path.join(d, f) for f in files]


def run(exe, files, kenv=None):
    env = dict(os.environ)
    if kenv:
        env["K"] = kenv
    t0 = time.monotonic()
    r = subprocess.run([os.path.join(ROOT, exe)] + files, capture_output=True, text=True, env=env)
    dt = time.monotonic() - t0
    rows = {}
    for line in r.stdout.splitlines():
        p = line.split()
        if len(p) != 7:
            continue
        rows[os.path.abspath(p[0])] = dict(w=float(p[1]), kpi=float(p[2]), T=int(p[3]), digs=int(p[4]),
                                           cpu=float(p[5]), ok=int(p[6]))
    return rows, dt, r.stderr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sol", action="append", required=True, help="файл решения (можно несколько)")
    ap.add_argument("--sets", nargs="*", help="папки с тестами (по умолчанию open + build/gen/*)")
    ap.add_argument("--K", default=None, help="ручки: '3=1.5,7=0'")
    ap.add_argument("--holdout", action="store_true", help="использовать build/holdout вместо build/gen")
    args = ap.parse_args()

    if args.sets:
        sets = []
        for s in args.sets:
            ap_ = os.path.abspath(s)
            if any(f.endswith(".in") for f in os.listdir(ap_)):
                sets.append(ap_)
            else:   # папка-контейнер: каждый подкаталог — отдельный набор
                sets += [os.path.join(ap_, d) for d in sorted(os.listdir(ap_))
                         if os.path.isdir(os.path.join(ap_, d))]
    else:
        sets = [os.path.join(ROOT, "materials", "tests")]
        base = os.path.join(ROOT, "build", "holdout" if args.holdout else "gen")
        if os.path.isdir(base):
            sets += [os.path.join(base, d) for d in sorted(os.listdir(base))
                     if os.path.isdir(os.path.join(base, d))]
    set_files = [(("open" if os.path.abspath(s) == os.path.join(ROOT, "materials", "tests")
                   else os.path.basename(os.path.abspath(s))), tests_in(s))
                 for s in sets if os.path.isdir(s)]
    all_files = [f for _, fs in set_files for f in fs]
    if not all_files:
        sys.exit("no tests found — run scripts/gen_tests.py first")

    results = {}
    for i, sol in enumerate(args.sol):
        out = f"build/sim_{i}.exe"
        build(sol, out)
        rows, dt, err = run(out, all_files, args.K)
        results[sol] = (rows, dt, err, out)
        if err.strip():
            print(f"{sol}: stderr:\n{err[:2000]}", file=sys.stderr)

    # таблица: строка = набор, колонки = решения
    labels = [os.path.basename(s) for s in args.sol]
    w = max(12, max(len(l) for l in labels) + 2)
    print(f"{'set':<22}{'n':>4} " + "".join(f"{l:>{w}}" for l in labels))
    print("-" * (26 + w * len(labels)))
    for name, files in set_files:
        line = f"{name:<22}{len(files):>4} "
        for sol in args.sol:
            rows = results[sol][0]
            vals = [rows[os.path.abspath(f)]["kpi"] for f in files if os.path.abspath(f) in rows]
            ok = all(rows.get(os.path.abspath(f), {}).get("ok", 0) for f in files)
            mean = sum(vals) / len(vals) if vals else 0.0
            line += f"{mean:>{w}.5f}" + ("" if ok else "!")
        print(line)
    # open: взвешенная сумма
    for name, files in set_files:
        if name != "open":
            continue
        line = f"{'open Sw':<22}{len(files):>4} "
        for sol in args.sol:
            rows = results[sol][0]
            tot = sum(rows[os.path.abspath(f)]["w"] * rows[os.path.abspath(f)]["kpi"]
                      for f in files if os.path.abspath(f) in rows)
            line += f"{tot:>{w}.4f}"
        print(line)
    print()
    for sol in args.sol:
        rows, dt, err, out = results[sol]
        cpu = sum(r["cpu"] for r in rows.values())
        digs = sum(r["digs"] for r in rows.values())
        # сводная метрика: open Σw + средний KPI по сгенерированным наборам (не open)
        opensum = sum(rows[os.path.abspath(f)]["w"] * rows[os.path.abspath(f)]["kpi"]
                      for name, files in set_files if name == "open" for f in files
                      if os.path.abspath(f) in rows)
        gen = [rows[os.path.abspath(f)]["kpi"] for name, files in set_files if name != "open"
               for f in files if os.path.abspath(f) in rows]
        genmean = sum(gen) / len(gen) if gen else 0.0
        print(f"{sol}: {len(rows)} tests, {dt:.1f}s wall, sum CPU {cpu:.2f}s, digs {digs}, "
              f"open {opensum:.4f}, gen-mean {genmean:.5f}")


if __name__ == "__main__":
    main()
