#!/usr/bin/env python3
"""Жадный подбор портфеля политик по CSV из `sim --pol all --pool <pool> --csv <out>`.

Метрика портфеля S — «oracle»: на каждом тесте берётся лучшая политика из S (верхняя оценка для выбора по истории),
затем среднее по классам с весами CLASS_W. Это не то, что даст реальный выбор, но показывает, какие политики дополняют друг друга.

  python scripts/portfolio_opt.py build/pool1.csv [--k 10]
"""
import argparse
import csv
import collections
import re

CLASS_W = {"open": 40.0, "twophase": 10.0, "alternating": 8.0, "sorted2": 3.0, "sortedbig": 6.0, "invert": 2.0, "periodic": 2.0, "phases": 4.0, "tiny": 1.0, "sawtooth": 2.0, "bimodal": 2.0, "fifo": 2.0, "lifo": 1.0, "desc": 2.0, "inv": 3.0, "scr": 3.0, "shift": 3.0, "half": 3.0, "dbl": 3.0}


def klass(label):
    if label.startswith("materials"):
        return "open"
    return label.split("/")[-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--k", type=int, default=10)
    ap.add_argument("--start", default="")
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.csv)))
    pols = sorted({r["pol"] for r in rows}, key=lambda p: [r["pol"] for r in rows].index(p))
    tests = {}
    for r in rows:
        key = (r["set"], r["test"])
        tests.setdefault(key, {"w": float(r["w"]), "cls": klass(r["set"]), "k": {}})["k"][r["pol"]] = float(r["kpi"])
    keys = list(tests)
    pidx = {p: i for i, p in enumerate(pols)}
    # матрица (чистый Python)
    cls = [tests[k]["cls"] for k in keys]
    w = [tests[k]["w"] for k in keys]
    M = [[tests[k]["k"][p] for p in pols] for k in keys]   # M[test][pol]
    classes = sorted(set(cls))
    idxs = {c: [i for i, x in enumerate(cls) if x == c] for c in classes}

    def per_class(vec):
        out = {}
        for c in classes:
            ii = idxs[c]
            out[c] = sum(vec[i] * w[i] for i in ii) if c == "open" else sum(vec[i] for i in ii) / len(ii)
        return out

    def score(vec):
        pc = per_class(vec)
        tot = 0.0
        for c in classes:
            v = pc[c] / sum(w[i] for i in idxs[c]) if c == "open" else pc[c]
            tot += CLASS_W.get(c, 1) * v
        return tot

    col = lambda i: [M[t][i] for t in range(len(keys))]
    singles = sorted(((score(col(i)), p) for p, i in pidx.items()), reverse=True)
    print("лучшие одиночные политики:")
    for s_, p in singles[:8]:
        print(f"  {s_:8.3f}  {p}")
    chosen = [p for p in a.start.split(",") if p]
    best = [0.0] * len(keys)
    for p in chosen:
        best = [max(x, y) for x, y in zip(best, col(pidx[p]))]
    cur = score(best) if chosen else 0.0
    print("\nжадный набор (oracle):")
    for step in range(a.k):
        gains = []
        for p, i in pidx.items():
            if p in chosen:
                continue
            gains.append((score([max(x, y) for x, y in zip(best, col(i))]), p))
        s_, p = max(gains)
        chosen.append(p)
        best = [max(x, y) for x, y in zip(best, col(pidx[p]))]
        pc = per_class(best)
        print(f"  +{p:28s} total={s_:8.3f} (+{s_ - cur:6.3f})  " + " ".join(f"{c}={pc[c]:.4f}" if c != "open" else f"open={pc[c]:.3f}" for c in classes))
        cur = s_


if __name__ == "__main__":
    main()
