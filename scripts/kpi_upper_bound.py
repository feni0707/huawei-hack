"""Тривиальная верхняя оценка KPI: живые в конце / 20·ceil(пик живых / 20).

  python3 scripts/kpi_upper_bound.py materials/tests            # по тестам + взвешенная сумма
  python3 scripts/kpi_upper_bound.py build/general/*            # несколько папок: среднее UB по папке
  python3 scripts/kpi_upper_bound.py --csv build/ub.csv DIR...  # путь,UB для bench/анализа разрыва
"""
import sys, os, math, glob, re


def ub_of(path):
    t = open(path).read().split()
    w, n = float(t[0]), int(t[1]); cmds = list(map(int, t[2:2 + n]))
    alive = {}; pid = 0; cur = 0; peak = 0
    for c in cmds:
        if c > 0:
            pid += 1; alive[pid] = c; cur += c; peak = max(peak, cur)
        else:
            cur -= alive.pop(-c)
    return w, cur / (20 * math.ceil(peak / 20)) if peak else 0.0


def key(p):
    m = re.search(r"(\d+)\.in$", p)
    return int(m.group(1)) if m else 0


args = sys.argv[1:]
csv = None
if args[:1] == ["--csv"]:
    csv, args = args[1], args[2:]
rows = []
for d in args:
    files = sorted(glob.glob(os.path.join(d, "*.in")), key=key)
    res = [(f, *ub_of(f)) for f in files]
    rows += res
    if len(args) == 1:
        tot = sum(w * b for _, w, b in res); totw = sum(w for _, w, _ in res)
        for f, w, b in res:
            print(f"test {os.path.basename(f):>7} w={w:<4} UB_KPI={b:.4f}  weighted={b * w:.3f}")
        print(f"sum UB={tot:.3f} of max {totw:g}")
    else:
        print(f"{os.path.basename(d.rstrip('/')):14s} n={len(res):3d} avg UB={sum(b for *_, b in res) / max(1, len(res)):.4f}")
if csv:
    with open(csv, "w") as f:
        for p, _, b in rows:
            f.write(f"{os.path.relpath(p)},{b:.6f}\n")
