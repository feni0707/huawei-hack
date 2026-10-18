#!/usr/bin/env python3
"""Бенчмарк решений: прогон нескольких бинарников на наборах тестов параллельно, без runner.py.

Каждый прогон — build/interactor (собирается scripts/build_interactor.sh) + решение, соединённые пайпами,
в отдельной временной папке (интерактор пишет score.txt в cwd). CPU решения меряется точно через wait4.
TL = CPU решения > --tl (5 с, как на NSUts) или превышен --timeout по wall-clock (зависание / нет flush).
WA и TL считаются как KPI = 0 (так же, как на проверке).

  python3 scripts/bench.py build/naive build/greedy            # таблица: open Total + avg KPI по классам
  python3 scripts/bench.py solutions/greedy.cpp                 # .cpp собирается в build/bin/<имя> (-O2)
  python3 scripts/bench.py --diff build/naive build/greedy      # потестовое сравнение двух решений
  python3 scripts/bench.py build/greedy --tests build/gen/random_bursts -j 4

Наборы по умолчанию: materials/tests (open) и каждая папка build/gen/* (см. scripts/gen_tests.py).
Для набора с весами ≠ 1 (open) выводится Σ w·KPI, иначе средний KPI. Сырые результаты: build/bench/last.csv.
CPU при параллельном прогоне завышается (конкуренция за ядра, пайпы); для сравнения CPU (тай-брейк) — -j 1.
"""
import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OPEN_DIR = os.path.join(ROOT, "materials", "tests")
GEN_DIR = os.path.join(ROOT, "build", "gen")
SCORE_RE = re.compile(r"Score = (\d+)/(\d+) \* weight = ([0-9.]+) x")

USE_COLOR = sys.stdout.isatty() and os.environ.get("NO_COLOR") is None


def color(s, code):
    return f"\033[{code}m{s}\033[0m" if USE_COLOR else s


def red(s):
    return color(s, "1;31")


def green(s):
    return color(s, "32")


# ---------------------------------------------------------------- подготовка

def build_interactor():
    out = subprocess.run([os.path.join(ROOT, "scripts", "build_interactor.sh")],
                         capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit("failed to build interactor:\n" + out.stderr)
    return out.stdout.strip()


def resolve_solution(path):
    """Бинарник как есть; .cpp — собрать в build/bin/<stem> (если исходник новее)."""
    if path.endswith(".cpp"):
        stem = os.path.splitext(os.path.basename(path))[0]
        exe = os.path.join(ROOT, "build", "bin", stem)
        os.makedirs(os.path.dirname(exe), exist_ok=True)
        if not os.path.exists(exe) or os.path.getmtime(path) > os.path.getmtime(exe):
            cxx = os.environ.get("CXX", "c++")
            print(f"compiling {path} -> {os.path.relpath(exe, ROOT)}", file=sys.stderr)
            r = subprocess.run([cxx, "-std=c++17", "-O2", "-o", exe, path])
            if r.returncode != 0:
                sys.exit(f"compilation of {path} failed")
        return os.path.abspath(exe)
    if not (os.path.isfile(path) and os.access(path, os.X_OK)):
        sys.exit(f"{path}: not an executable file (pass a binary or a .cpp)")
    return os.path.abspath(path)


def test_sort_key(name):
    stem = os.path.splitext(name)[0]
    return (0, int(stem), "") if stem.isdigit() else (1, 0, stem)


def collect_sets(dirs):
    """[(label, [test paths])]; label 'open' для materials/tests, иначе имя папки."""
    sets = []
    for d in dirs:
        d = os.path.abspath(d)
        if not os.path.isdir(d):
            sys.exit(f"{d}: not a directory")
        files = sorted((f for f in os.listdir(d) if f.endswith(".in")), key=test_sort_key)
        if not files:
            print(f"warning: no *.in in {d}", file=sys.stderr)
            continue
        label = "open" if os.path.samefile(d, OPEN_DIR) else os.path.basename(d)
        if label in (l for l, _ in sets):   # одноимённые папки из разных мест
            label = os.path.relpath(d, ROOT)
        sets.append((label, [os.path.join(d, f) for f in files]))
    return sets


def read_weight(path):
    with open(path) as f:
        return float(f.readline().split()[0])


# ---------------------------------------------------------------- один прогон

def run_one(interactor, sol, test, timeout, tl):
    """-> dict(status, kpi, occupied, trench, cpu, wall, msg)."""
    res = dict(status="OK", kpi=0.0, occupied=0, trench=0, cpu=0.0, wall=0.0, msg="")
    with tempfile.TemporaryDirectory(prefix="bench_") as work:
        i2s_r, i2s_w = os.pipe()   # интерактор -> решение
        s2i_r, s2i_w = os.pipe()   # решение -> интерактор
        t0 = time.monotonic()
        inter = subprocess.Popen([interactor, test, os.path.join(work, "tout.txt")],
                                 stdin=s2i_r, stdout=i2s_w, stderr=subprocess.PIPE, cwd=work)
        solp = subprocess.Popen([sol], stdin=i2s_r, stdout=s2i_w, stderr=subprocess.DEVNULL, cwd=work,
                                start_new_session=True)   # своя группа: на таймауте убиваем и потомков
        for fd in (i2s_r, i2s_w, s2i_r, s2i_w):
            os.close(fd)

        timed_out = threading.Event()

        def kill():
            timed_out.set()
            try:
                os.killpg(solp.pid, 9)
            except OSError:
                pass
            try:
                inter.kill()
            except OSError:
                pass

        timer = threading.Timer(timeout, kill)
        timer.start()
        try:
            err = inter.communicate()[1].decode(errors="replace")
            _, status, ru = os.wait4(solp.pid, 0)
            solp.returncode = os.waitstatus_to_exitcode(status)   # чтобы Popen не ждал повторно
        finally:
            timer.cancel()
        res["wall"] = time.monotonic() - t0
        res["cpu"] = ru.ru_utime + ru.ru_stime

    m = SCORE_RE.search(err)
    if timed_out.is_set():
        res["status"], res["msg"] = "TL", f"wall-clock timeout {timeout:g}s"
    elif m:
        res["occupied"], res["trench"] = int(m.group(1)), int(m.group(2))
        res["kpi"] = res["occupied"] / res["trench"]
        if res["cpu"] > tl:
            res["status"], res["msg"] = "TL", f"CPU {res['cpu']:.2f}s > {tl:g}s"
    else:
        lines = [l for l in err.strip().splitlines() if l.strip()]
        res["status"] = "WA"
        res["msg"] = lines[-1] if lines else f"interactor exit {inter.returncode}, no output"
        if solp.returncode not in (0, None) and solp.returncode < 0:
            res["msg"] += f" (solution killed by signal {-solp.returncode})"
    if res["status"] != "OK":
        res["kpi"] = 0.0
    return res


# ---------------------------------------------------------------- отчёт

def fmt_table(header, rows, align=None):
    plain = lambda s: re.sub(r"\033\[[0-9;]*m", "", s)
    widths = [max(len(plain(str(r[i]))) for r in [header] + rows) for i in range(len(header))]
    align = align or ["<"] + [">"] * (len(header) - 1)
    out = []
    for k, r in enumerate([header] + rows):
        cells = []
        for i, c in enumerate(r):
            pad = widths[i] - len(plain(str(c)))
            cells.append(str(c) + " " * pad if align[i] == "<" else " " * pad + str(c))
        out.append("  ".join(cells))
        if k == 0:
            out.append("  ".join("-" * w for w in widths))
    return "\n".join(out)


def set_metric(label, tests, results, sol):
    """(значение, подпись): Σ w·KPI для набора с весами ≠ 1, иначе средний KPI."""
    ws = [read_weight(t) for t in tests]
    kpis = [results[(sol, t)]["kpi"] for t in tests]
    if any(w != 1 for w in ws):
        return sum(w * k for w, k in zip(ws, kpis)), f"{label} Σw (/{sum(ws):g})"
    return sum(kpis) / len(kpis), f"{label} avgKPI"


def summary(sols, labels, sets, results):
    header, rows = ["solution"], []
    for li, sol in enumerate(sols):
        row = [labels[li]]
        for label, tests in sets:
            v, h = set_metric(label, tests, results, sol)
            if li == 0:
                header.append(h)
            bad = any(results[(sol, t)]["status"] != "OK" for t in tests)
            row.append(red(f"{v:.4f}") if bad else f"{v:.4f}")
        all_res = [results[(sol, t)] for _, tests in sets for t in tests]
        wa = sum(r["status"] == "WA" for r in all_res)
        tl = sum(r["status"] == "TL" for r in all_res)
        mx = max(all_res, key=lambda r: r["cpu"])
        row += [red(str(wa)) if wa else "0", red(str(tl)) if tl else "0",
                f"{mx['cpu']:.3f} ({mx['test']})", f"{sum(r['cpu'] for r in all_res):.2f}"]
        rows.append(row)
    header += ["WA", "TL", "max CPU s (test)", "Σ CPU s"]
    return fmt_table(header, rows)


def diff_report(a, b, la, lb, sets, results):
    rows, wins, losses = [], 0, 0
    for label, tests in sets:
        for t in tests:
            ra, rb = results[(a, t)], results[(b, t)]
            d = rb["kpi"] - ra["kpi"]
            if d > 1e-9:
                wins += 1
                ds = green(f"{d:+.4f}")
            elif d < -1e-9:
                losses += 1
                ds = red(f"{d:+.4f}")
            else:
                ds = f"{d:+.4f}"
            sa = f"{ra['kpi']:.4f}" if ra["status"] == "OK" else red(ra["status"])
            sb = f"{rb['kpi']:.4f}" if rb["status"] == "OK" else red(rb["status"])
            rows.append([f"{label}/{os.path.basename(t)}", f"{read_weight(t):g}", sa, sb, ds,
                         f"{ra['cpu']:.3f}", f"{rb['cpu']:.3f}"])
    header = ["test", "w", f"KPI {la}", f"KPI {lb}", "Δ (B−A)", f"CPU {la}", f"CPU {lb}"]
    return fmt_table(header, rows) + f"\n\n{lb} vs {la}: better on {wins}, worse on {losses}, " \
                                    f"equal on {len(rows) - wins - losses} tests"


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("solutions", nargs="*", help="бинарники решений или .cpp")
    ap.add_argument("--diff", nargs=2, metavar=("A", "B"), help="потестовое сравнение двух решений")
    ap.add_argument("--tests", nargs="+", metavar="DIR",
                    help="папки с *.in (по умолчанию materials/tests и build/gen/*)")
    ap.add_argument("-j", "--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2),
                    help="параллельных прогонов (по умолчанию ядра/2); для точного CPU — -j 1")
    ap.add_argument("--timeout", type=float, default=30.0, help="wall-clock лимит на прогон, с (30)")
    ap.add_argument("--tl", type=float, default=5.0, help="лимит CPU решения, с (5, как на NSUts)")
    args = ap.parse_args()

    paths = list(args.diff) if args.diff else args.solutions
    if not paths:
        ap.error("give solution binaries or --diff A B")
    if args.diff and args.solutions:
        ap.error("--diff takes exactly two solutions; don't pass extra positional ones")

    dirs = args.tests
    if not dirs:
        dirs = [OPEN_DIR]
        if os.path.isdir(GEN_DIR):
            dirs += sorted(os.path.join(GEN_DIR, d) for d in os.listdir(GEN_DIR)
                           if os.path.isdir(os.path.join(GEN_DIR, d)))
        if len(dirs) == 1:
            print("note: build/gen is empty — run scripts/gen_tests.py for generated sets", file=sys.stderr)
    sets = collect_sets(dirs)

    interactor = build_interactor()
    sols = [resolve_solution(p) for p in paths]
    names = [os.path.basename(p) for p in sols]
    labels = names if len(set(names)) == len(names) else [os.path.relpath(p, ROOT) for p in sols]

    jobs = [(s, t) for s in sols for _, tests in sets for t in tests]
    # Тяжёлые тесты — первыми, чтобы не ждать их в хвосте.
    jobs.sort(key=lambda j: -os.path.getsize(j[1]))
    label_of = {t: f"{lab}/{os.path.basename(t)}" for lab, tests in sets for t in tests}
    results = {}
    t0 = time.monotonic()
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futs = {ex.submit(run_one, interactor, s, t, args.timeout, args.tl): (s, t) for s, t in jobs}
        for k, f in enumerate(as_completed(futs), 1):
            s, t = futs[f]
            r = f.result()
            r["test"] = label_of[t]
            results[(s, t)] = r
            if r["status"] != "OK":
                print(red(f"\n{r['status']}: {labels[sols.index(s)]} on {r['test']}: {r['msg']}"),
                      file=sys.stderr)
            if sys.stderr.isatty():
                print(f"\r[{k}/{len(jobs)}] {time.monotonic() - t0:.1f}s", end="", file=sys.stderr, flush=True)
    if sys.stderr.isatty():
        print(file=sys.stderr)

    os.makedirs(os.path.join(ROOT, "build", "bench"), exist_ok=True)
    with open(os.path.join(ROOT, "build", "bench", "last.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["solution", "set", "test", "weight", "status", "kpi", "occupied", "trench", "cpu", "wall", "msg"])
        for li, s in enumerate(sols):
            for lab, tests in sets:
                for t in tests:
                    r = results[(s, t)]
                    w.writerow([labels[li], lab, os.path.basename(t), read_weight(t), r["status"],
                                f"{r['kpi']:.6f}", r["occupied"], r["trench"], f"{r['cpu']:.4f}",
                                f"{r['wall']:.3f}", r["msg"]])

    if args.diff:
        print(diff_report(sols[0], sols[1], labels[0], labels[1], sets, results))
        print()
    print(summary(sols, labels, sets, results))
    bad = [(labels[sols.index(s)], r) for (s, t), r in results.items() if r["status"] != "OK"]
    if bad:
        print("\n" + red("FAILURES:"))
        for lab, r in sorted(bad, key=lambda x: (x[0], x[1]["test"])):
            print(red(f"  {r['status']}  {lab}  {r['test']}: {r['msg']}"))
    print(f"\n{len(jobs)} runs in {time.monotonic() - t0:.1f}s, -j {args.jobs}; raw: build/bench/last.csv")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
