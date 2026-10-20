#!/usr/bin/env python3
"""Windows-порт scripts/bench.py: прогон решений через официальный интерактор, параллельно.

Отличия от bench.py (POSIX-only: os.pipe+wait4+killpg):
  * CPU решения читается через WinAPI GetProcessTimes (опрос из потока, ctypes, без psutil);
  * таймаут убивает дерево процессов через taskkill /T /F;
  * .cpp собирается MSVC (cl.exe) по vcvars64.bat, кэш в build/bin/<stem>.exe (без пересборки, если exe свежее);
  * KPI берётся ИЗ score.txt, который пишет интерактор (а не из regex по stderr) — это то же число,
    что уходит в проверяющую систему, и заодно проверяет, что интерактор дошёл до конца;
  * stderr интерактора сохраняется при WA (причина отказа).

  python scripts/winbench.py build/bin/naive.exe
  python scripts/winbench.py solutions/foo.cpp                  # соберёт и прогонит
  python scripts/winbench.py build/bin/foo.exe --tests build/gen/two_phase_real -j 8
  python scripts/winbench.py --diff build/bin/a.exe build/bin/b.exe
  python scripts/winbench.py build/bin/foo.exe -j 1             # точный CPU

Наборы по умолчанию: materials/tests (open) + build/gen/*. Для набора с весами != 1 — Σ w·KPI.
Сырые результаты: build/bench/win_last.csv
"""
import argparse
import csv
import ctypes
import ctypes.wintypes as wt
import itertools
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OPEN_DIR = os.path.join(ROOT, "materials", "tests")
GEN_DIR = os.path.join(ROOT, "build", "gen")
INTERACTOR = os.path.join(ROOT, "materials", "interactor.exe")
BENCH_DIR = os.path.join(ROOT, "build", "bench")

USE_COLOR = sys.stdout.isatty() and os.environ.get("NO_COLOR") is None


def color(s, code):
    return f"\033[{code}m{s}\033[0m" if USE_COLOR else s


def red(s):
    return color(s, "1;31")


def green(s):
    return color(s, "32")


# ------------------------------------------------------------------ WinAPI CPU
_k32 = ctypes.WinDLL("kernel32", use_last_error=True)
_k32.OpenProcess.restype = wt.HANDLE
_k32.OpenProcess.argtypes = [wt.DWORD, wt.BOOL, wt.DWORD]
_k32.GetProcessTimes.argtypes = [wt.HANDLE] + [ctypes.POINTER(wt.FILETIME)] * 4
_k32.CloseHandle.argtypes = [wt.HANDLE]
PROCESS_QUERY_LIMITED_INFORMATION = 0x1000


class CpuSampler(threading.Thread):
    """Опрос GetProcessTimes решения, пока процесс жив; последнее значение сохраняется."""

    def __init__(self, pid, interval=0.01):
        super().__init__(daemon=True)
        self.h = _k32.OpenProcess(PROCESS_QUERY_LIMITED_INFORMATION, False, pid)
        self.interval = interval
        self.cpu = 0.0
        self.stop = threading.Event()

    def read(self):
        if not self.h:
            return 0.0
        c, e, k, u = (wt.FILETIME() for _ in range(4))
        if not _k32.GetProcessTimes(self.h, ctypes.byref(c), ctypes.byref(e),
                                    ctypes.byref(k), ctypes.byref(u)):
            return self.cpu
        val = lambda ft: (ft.dwHighDateTime << 32 | ft.dwLowDateTime) / 1e7
        return val(k) + val(u)

    def run(self):
        while not self.stop.is_set():
            self.cpu = max(self.cpu, self.read())
            self.stop.wait(self.interval)
        self.cpu = max(self.cpu, self.read())

    def finish(self):
        """Остановить опрос и дочитать финальные значения (handle жив, пока не закрыт)."""
        self.stop.set()
        self.join(timeout=2)
        for _ in range(50):                    # процесс уже мёртв: значения окончательные
            v = self.read()
            if v == self.cpu:
                break
            self.cpu = max(self.cpu, v)
            time.sleep(0.002)
        if self.h:
            _k32.CloseHandle(self.h)
            self.h = None
        return self.cpu


# ------------------------------------------------------------------ сборка
def find_vcvars():
    vswhere = os.path.join(os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)"),
                           "Microsoft Visual Studio", "Installer", "vswhere.exe")
    if os.path.exists(vswhere):
        try:
            out = subprocess.run([vswhere, "-latest", "-property", "installationPath"],
                                 capture_output=True, text=True, timeout=30).stdout.strip()
            if out:
                p = os.path.join(out, "VC", "Auxiliary", "Build", "vcvars64.bat")
                if os.path.exists(p):
                    return p
        except Exception:
            pass
    for cand in (r"C:\Program Files\Microsoft Visual Studio\2022\Community",
                 r"C:\Program Files\Microsoft Visual Studio\2022\Professional",
                 r"C:\Program Files\Microsoft Visual Studio\2022\Enterprise"):
        p = os.path.join(cand, "VC", "Auxiliary", "Build", "vcvars64.bat")
        if os.path.exists(p):
            return p
    sys.exit("vcvars64.bat not found (нужен MSVC 2022)")


VCVARS = None


def build_cpp(path):
    """cl /O2 в build/bin/<stem>.exe; пересборка, если .cpp новее."""
    global VCVARS
    stem = os.path.splitext(os.path.basename(path))[0]
    exe = os.path.join(ROOT, "build", "bin", stem + ".exe")
    os.makedirs(os.path.dirname(exe), exist_ok=True)
    if os.path.exists(exe) and os.path.getmtime(path) <= os.path.getmtime(exe):
        return exe
    if VCVARS is None:
        VCVARS = find_vcvars()
    cmd = (f'call "{VCVARS}" >nul 2>nul && cd /d "{ROOT}" && '
           f'cl /nologo /O2 /EHsc /std:c++17 /Fo:build\\bin\\ /Fe:build\\bin\\{stem}.exe "{path}"')
    print(f"compiling {os.path.relpath(path, ROOT)} -> build/bin/{stem}.exe", file=sys.stderr)
    r = subprocess.run(["cmd", "/c", cmd], capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(exe):
        sys.exit(f"compilation failed:\n{r.stdout}\n{r.stderr}")
    return exe


def resolve_solution(path):
    if path.endswith(".cpp"):
        return build_cpp(os.path.abspath(path))
    if not os.path.isfile(path):
        sys.exit(f"{path}: not a file")
    return os.path.abspath(path)


# ------------------------------------------------------------------ тесты
def test_sort_key(name):
    stem = os.path.splitext(name)[0]
    return (0, int(stem), "") if stem.isdigit() else (1, 0, stem)


def collect_sets(dirs):
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
        if label in (l for l, _ in sets):
            label = os.path.relpath(d, ROOT)
        sets.append((label, [os.path.join(d, f) for f in files]))
    return sets


WEIGHT_CACHE = {}


def read_weight(path):
    if path not in WEIGHT_CACHE:
        with open(path) as f:
            WEIGHT_CACHE[path] = float(f.readline().split()[0])
    return WEIGHT_CACHE[path]


# ------------------------------------------------------------------ прогон
RUN_DIR = os.path.join(ROOT, "build", "run")
_work_ids = itertools.count(1)
# ВАЖНО: рабочие папки только внутри репозитория. Официальный interactor.exe отказывается создавать
# test-output-file в %TEMP% ("Can not write to the test-output-file") — проверено на этой машине.
_WORK_LOCK = threading.Lock()


def make_work():
    with _WORK_LOCK:
        d = os.path.join(RUN_DIR, f"wb_{os.getpid()}_{next(_work_ids)}")
    os.makedirs(d, exist_ok=True)
    return d


def run_one(sol, test, timeout, tl):
    res = dict(status="OK", kpi=0.0, occupied=0, trench=0, cpu=0.0, wall=0.0, msg="")
    work = make_work()
    score_path = os.path.join(work, "score.txt")
    i2s_r, i2s_w = os.pipe()          # интерактор -> решение
    s2i_r, s2i_w = os.pipe()          # решение -> интерактор
    for fd in (i2s_r, i2s_w, s2i_r, s2i_w):
        os.set_inheritable(fd, True)
    t0 = time.monotonic()
    inter = solp = None
    try:
        inter = subprocess.Popen([INTERACTOR, test, os.path.join(work, "tout.txt")],
                                 stdin=s2i_r, stdout=i2s_w, stderr=subprocess.PIPE, cwd=work)
        solp = subprocess.Popen([sol], stdin=i2s_r, stdout=s2i_w,
                                stderr=subprocess.DEVNULL, cwd=work)
        for fd in (i2s_r, i2s_w, s2i_r, s2i_w):
            os.close(fd)
        i2s_r = i2s_w = s2i_r = s2i_w = -1
        sampler = CpuSampler(solp.pid)
        sampler.start()
        timed_out = False
        try:
            err = inter.communicate(timeout=timeout)[1].decode(errors="replace")
        except subprocess.TimeoutExpired:
            timed_out = True
            err = ""
            subprocess.run(["taskkill", "/T", "/F", "/PID", str(inter.pid)],
                           capture_output=True)
            try:
                subprocess.run(["taskkill", "/T", "/F", "/PID", str(solp.pid)],
                               capture_output=True)
            except Exception:
                pass
        try:
            solp.wait(timeout=10)
        except subprocess.TimeoutExpired:
            subprocess.run(["taskkill", "/T", "/F", "/PID", str(solp.pid)], capture_output=True)
        res["cpu"] = sampler.finish()
        res["wall"] = time.monotonic() - t0
    finally:
        for fd in (i2s_r, i2s_w, s2i_r, s2i_w):
            if fd >= 0:
                try:
                    os.close(fd)
                except OSError:
                    pass

    # KPI — из score.txt (его пишет интерактор: KPI*weight); occupied/total — из stderr
    m = re.search(r"Score = (\d+)/(\d+) \* weight", err)
    if m:
        res["occupied"], res["trench"] = int(m.group(1)), int(m.group(2))
        res["kpi"] = res["occupied"] / res["trench"] if res["trench"] else 0.0
    if timed_out:
        res["status"], res["msg"] = "TL", f"wall timeout {timeout:g}s"
    elif not m:
        lines = [l for l in err.strip().splitlines() if l.strip()]
        res["status"] = "WA"
        res["msg"] = lines[-1] if lines else f"interactor exit {inter.returncode}, no score"
    elif res["cpu"] > tl:
        res["status"], res["msg"] = "TL", f"CPU {res['cpu']:.2f}s > {tl:g}s"
    if res["status"] != "OK":
        res["kpi"] = 0.0
    shutil.rmtree(work, ignore_errors=True)
    return res


# ------------------------------------------------------------------ отчёт
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
    ws = [read_weight(t) for t in tests]
    kpis = [results[(sol, t)]["kpi"] for t in tests]
    if any(w != 1 for w in ws):
        return sum(w * k for w, k in zip(ws, kpis)), f"{label} Sw (/{sum(ws):g})"
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
    header += ["WA", "TL", "max CPU s (test)", "Sum CPU s"]
    return fmt_table(header, rows)


def main():
    global INTERACTOR
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("solutions", nargs="*")
    ap.add_argument("--diff", nargs=2, metavar=("A", "B"))
    ap.add_argument("--tests", nargs="+", metavar="DIR")
    ap.add_argument("-j", "--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--timeout", type=float, default=30.0)
    ap.add_argument("--tl", type=float, default=5.0)
    ap.add_argument("--csv", default=os.path.join(BENCH_DIR, "win_last.csv"))
    ap.add_argument("--interactor", default=INTERACTOR, help="путь к interactor.exe")
    args = ap.parse_args()

    paths = list(args.diff) if args.diff else args.solutions
    if not paths:
        ap.error("give solution binaries or --diff A B")

    dirs = args.tests
    if not dirs:
        dirs = [OPEN_DIR]
        if os.path.isdir(GEN_DIR):
            dirs += sorted(os.path.join(GEN_DIR, d) for d in os.listdir(GEN_DIR)
                           if os.path.isdir(os.path.join(GEN_DIR, d)))
    sets = collect_sets(dirs)
    INTERACTOR = os.path.abspath(args.interactor)
    if not os.path.exists(INTERACTOR):
        sys.exit(f"no interactor at {INTERACTOR}")

    sols = [resolve_solution(p) for p in paths]
    names = [os.path.basename(p) for p in sols]
    labels = names if len(set(names)) == len(names) else [os.path.relpath(p, ROOT) for p in sols]

    jobs = [(s, t) for s in sols for _, tests in sets for t in tests]
    jobs.sort(key=lambda j: -os.path.getsize(j[1]))
    label_of = {t: f"{lab}/{os.path.basename(t)}" for lab, tests in sets for t in tests}
    results = {}
    t0 = time.monotonic()
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futs = {ex.submit(run_one, s, t, args.timeout, args.tl): (s, t) for s, t in jobs}
        for k, f in enumerate(as_completed(futs), 1):
            s, t = futs[f]
            r = f.result()
            r["test"] = label_of[t]
            results[(s, t)] = r
            if r["status"] != "OK":
                print(red(f"\n{r['status']}: {labels[sols.index(s)]} on {r['test']}: {r['msg']}"),
                      file=sys.stderr)

    os.makedirs(BENCH_DIR, exist_ok=True)
    with open(args.csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["solution", "set", "test", "weight", "status", "kpi", "occupied",
                    "trench", "cpu", "wall", "msg"])
        for li, s in enumerate(sols):
            for lab, tests in sets:
                for t in tests:
                    r = results[(s, t)]
                    w.writerow([labels[li], lab, os.path.basename(t), read_weight(t), r["status"],
                                f"{r['kpi']:.6f}", r["occupied"], r["trench"], f"{r['cpu']:.4f}",
                                f"{r['wall']:.3f}", r["msg"]])

    if args.diff:
        a, b = sols
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
                rows.append([f"{label}/{os.path.basename(t)}", f"{read_weight(t):g}",
                             f"{ra['kpi']:.4f}" if ra["status"] == "OK" else red(ra["status"]),
                             f"{rb['kpi']:.4f}" if rb["status"] == "OK" else red(rb["status"]),
                             ds, f"{ra['cpu']:.3f}", f"{rb['cpu']:.3f}"])
        print(fmt_table(["test", "w", f"KPI {labels[0]}", f"KPI {labels[1]}", "D (B-A)",
                         f"CPU {labels[0]}", f"CPU {labels[1]}"], rows))
        print(f"\n{labels[1]} vs {labels[0]}: better {wins}, worse {losses}, "
              f"equal {len(rows) - wins - losses}\n")
    print(summary(sols, labels, sets, results))
    bad = [(labels[sols.index(s)], r) for (s, t), r in results.items() if r["status"] != "OK"]
    if bad:
        print("\n" + red("FAILURES:"))
        for lab, r in sorted(bad, key=lambda x: (x[0], x[1]["test"])):
            print(red(f"  {r['status']}  {lab}  {r['test']}: {r['msg']}"))
    print(f"\n{len(jobs)} runs in {time.monotonic() - t0:.1f}s, -j {args.jobs}; raw: "
          f"{os.path.relpath(args.csv, ROOT)}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
