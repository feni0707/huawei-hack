#!/usr/bin/env python3
"""Прогон готового .exe решения через настоящий интерактор (materials/interactor.exe) на Windows, параллельно.

Как runner.py, но без psutil: CPU решения берётся из GetProcessTimes. Корректность (WA) проверяет сам интерактор.

  python scripts/win_bench.py build/bin/portfolio.exe materials/tests build/gen2/twophase
  python scripts/win_bench.py build/bin/portfolio.exe materials/tests -j 4 -v

Для папок с весами ≠ 1 (open) печатается Σ w·KPI, иначе средний KPI. WA/Timeout = KPI 0.
"""
import argparse
import ctypes
import ctypes.wintypes as wt
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INTERACTOR = os.path.join(ROOT, "materials", "interactor.exe")


def cpu_seconds(handle):
    c, e, k, u = wt.FILETIME(), wt.FILETIME(), wt.FILETIME(), wt.FILETIME()
    ok = ctypes.windll.kernel32.GetProcessTimes(wt.HANDLE(int(handle)), ctypes.byref(c), ctypes.byref(e), ctypes.byref(k), ctypes.byref(u))
    if not ok:
        return -1.0
    t = lambda f: ((f.dwHighDateTime << 32) | f.dwLowDateTime) / 1e7
    return t(k) + t(u)


def run_one(exe, test, timeout):
    with tempfile.TemporaryDirectory() as d:
        score, ans, rep = (os.path.join(d, x) for x in ("score.txt", "answer.txt", "report.txt"))
        open(ans, "w").close()
        try:
            ip = subprocess.Popen([INTERACTOR, test, score, ans, rep], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE, cwd=d)
            sp = subprocess.Popen([exe], stdin=ip.stdout, stdout=ip.stdin, cwd=d)
            try:
                sp.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                sp.kill(); ip.kill()
                return 0.0, timeout, "TIMEOUT"
            cpu = cpu_seconds(sp._handle)
            try:
                ip.wait(timeout=5)
            except subprocess.TimeoutExpired:
                ip.kill()
            err = ip.stderr.read().decode("utf-8", "replace") if ip.stderr else ""
            ip.stdin.close() if ip.stdin and not ip.stdin.closed else None
        except Exception as ex:  # noqa: BLE001
            return 0.0, 0.0, f"ERR {ex}"
        try:
            sc = float(open(score).readline())
        except Exception:  # noqa: BLE001
            sc = 0.0
        m = re.search(r"Score = (\d+)/(\d+)", err)
        msg = err.strip().splitlines()[-1] if err.strip() else ""
        kpi = (int(m.group(1)) / int(m.group(2))) if m else 0.0
        return kpi if sc > 0 else 0.0, cpu, msg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("exe")
    ap.add_argument("paths", nargs="+")
    ap.add_argument("-j", type=int, default=max(1, (os.cpu_count() or 4) // 2))
    ap.add_argument("-v", action="store_true")
    ap.add_argument("--timeout", type=float, default=30)
    a = ap.parse_args()
    exe = os.path.abspath(a.exe)
    for p in a.paths:
        files = [os.path.abspath(p)] if p.endswith(".in") else sorted((os.path.abspath(os.path.join(p, f)) for f in os.listdir(p) if f.endswith(".in")),
                                                    key=lambda x: (len(os.path.basename(x)), x))
        weights = []
        for f in files:
            weights.append(float(open(f).read(40).split()[0]))
        with ThreadPoolExecutor(a.j) as ex:
            res = list(ex.map(lambda f: run_one(exe, f, a.timeout), files))
        weighted = any(w != 1 for w in weights)
        tot = sum(r[0] * w for r, w in zip(res, weights))
        bad = [(os.path.basename(f), r[2]) for f, r in zip(files, res) if r[0] == 0.0]
        name = os.path.basename(os.path.normpath(p))
        maxcpu = max(r[1] for r in res)
        if a.v:
            for f, r, w in zip(files, res, weights):
                print(f"  {os.path.basename(f):8s} w={w:<4g} kpi={r[0]:.4f} cpu={r[1]:.2f}s {r[2] if r[0] == 0 else ''}")
        print(f"{name:16s} {'sum' if weighted else 'avg'}={tot if weighted else tot / len(files):.4f}  maxcpu={maxcpu:.2f}s  failures={len(bad)} {bad[:3] if bad else ''}")


if __name__ == "__main__":
    main()
