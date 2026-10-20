#!/usr/bin/env python3
"""Разрыв до тривиальной верхней оценки по семействам (быстрый стенд scripts/simc.cpp, без интерактора).

  python3 scripts/gapeval.py build/simc_sf build/general/*                 # одно решение
  python3 scripts/gapeval.py build/simc_sf,build/simc_new build/gtrain/*   # сравнение: Δ KPI B−A по семействам

Для каждого семейства: средний KPI, средний UB (живые в конце / 20·ceil(пик/20)), эффективность Σkpi/Σub.
Для пары — средняя разница, ст. ошибка, победы/поражения, худший тест. Бинарники: «путь» или «путь|K-ручки».
J — потоков (по умолчанию 8). Тесты — папки с *.in; семейство = имя папки.
"""
import glob, math, os, re, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(__file__))
J = int(os.environ.get("J", "8"))


def ub_of(path):
    t = open(path).read().split()
    n = int(t[1]); alive = {}; pid = cur = peak = 0
    for c in map(int, t[2:2 + n]):
        if c > 0:
            pid += 1; alive[pid] = c; cur += c; peak = max(peak, cur)
        else:
            cur -= alive.pop(-c)
    return cur / (20 * math.ceil(peak / 20)) if peak else 0.0


def run(spec, files):
    b, _, k = spec.partition("|")
    out = subprocess.run([b] + files, env=dict(os.environ, K=k), capture_output=True, text=True).stdout
    r = {}
    for m in re.finditer(r"^(\S+)\s+w=\S+\s+kpi=(\S+).*$", out, re.M):
        r[m.group(1)] = (float(m.group(2)), "WA" in m.group(0))
    return r


def runall(spec, files):
    chunks = [files[i::J * 2] for i in range(J * 2)]
    res = {}
    with ThreadPoolExecutor(J) as ex:
        for r in ex.map(lambda c: run(spec, c) if c else {}, chunks):
            res.update(r)
    return res


def main():
    specs = sys.argv[1].split(",")
    dirs = sys.argv[2:]
    fam = {}
    for d in dirs:
        fs = sorted(glob.glob(os.path.join(d, "*.in")))
        if fs:
            fam[os.path.basename(d.rstrip("/"))] = fs
    files = [f for fs in fam.values() for f in fs]
    with ThreadPoolExecutor(J) as ex:
        ubs = dict(zip(files, ex.map(ub_of, files)))
    res = [runall(s, files) for s in specs]
    wa = [sum(1 for f in files if f not in r or r[f][1]) for r in res]
    hdr = f"{'family':14s} {'n':>3s} {'UB':>7s}" + "".join(f" {'KPI' + str(i):>8s} {'eff' + str(i):>7s}" for i in range(len(specs)))
    if len(specs) == 2:
        hdr += f" {'Δ(B−A)':>9s} {'±se':>7s} {'W/L':>7s}  worst"
    print(hdr)
    tot = [[0.0, 0.0] for _ in specs]; tub = 0.0; dsum = []
    for name, fs in fam.items():
        ub = sum(ubs[f] for f in fs)
        line = f"{name:14s} {len(fs):3d} {ub / len(fs):7.4f}"
        for i, r in enumerate(res):
            k = sum(r.get(f, (0, 1))[0] for f in fs)
            tot[i][0] += k
            line += f" {k / len(fs):8.4f} {k / ub if ub else 0:7.4f}"
        tub += ub
        if len(specs) == 2:
            d = [res[1].get(f, (0, 1))[0] - res[0].get(f, (0, 1))[0] for f in fs]
            dsum += d
            m = sum(d) / len(d)
            se = math.sqrt(sum((x - m) ** 2 for x in d) / max(1, len(d) - 1) / len(d))
            w = sum(x > 1e-9 for x in d); l = sum(x < -1e-9 for x in d)
            wi = min(range(len(d)), key=lambda i: d[i])
            line += f" {m:+9.5f} {se:7.5f} {w:3d}/{l:<3d}  {os.path.basename(fs[wi])} {d[wi]:+.4f}"
        print(line)
    n = len(files)
    line = f"{'ALL':14s} {n:3d} {tub / n:7.4f}" + "".join(f" {t[0] / n:8.4f} {t[0] / tub:7.4f}" for t in tot)
    if len(specs) == 2:
        m = sum(dsum) / n
        se = math.sqrt(sum((x - m) ** 2 for x in dsum) / max(1, n - 1) / n)
        line += f" {m:+9.5f} {se:7.5f} {sum(x > 1e-9 for x in dsum):3d}/{sum(x < -1e-9 for x in dsum):<3d}"
    print(line)
    if any(wa):
        print("WA/missing:", wa)


if __name__ == "__main__":
    main()
