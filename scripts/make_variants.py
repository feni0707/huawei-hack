#!/usr/bin/env python3
"""Варианты открытых тестов с другой «геометрией»: та же временная структура, другие длины труб.
Нужны, чтобы проверить, что решение учится на сценарии, а не подогнано под конкретные числа.
Результат: build/var/<map>/<N>.in (N = номер открытого теста)."""
import os, sys
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAPS = {
    "inv": lambda L: 41 - L,
    "scr": lambda L: (L * 7 - 1) % 40 + 1,
    "shift": lambda L: min(40, L + 4),
    "half": lambda L: max(1, L // 2),
    "dbl": lambda L: min(40, 2 * L),
}
TESTS = [5, 7, 8, 9, 10, 11, 12, 13, 14, 15]
for name, f in MAPS.items():
    d = os.path.join(root, "build", "var", name)
    os.makedirs(d, exist_ok=True)
    for t in TESTS:
        tok = open(os.path.join(root, "materials", "tests", f"{t}.in")).read().split()
        w, n = tok[0], int(tok[1])
        cmds = [f(int(x)) if int(x) > 0 else int(x) for x in tok[2:2 + n]]
        with open(os.path.join(d, f"{t}.in"), "w") as fo:
            fo.write(f"1 {n}\n" + "\n".join(map(str, cmds)) + "\n0\n")
    print(name, "->", d)
