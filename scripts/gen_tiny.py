#!/usr/bin/env python3
"""Генератор «крошечных» тестов-головоломок (как открытые 1, 2, 3, 6 и отчасти 4).

Смысл: на таких тестах KPI дискретен (T кратно 20, а суммарная длина мала), поэтому решение либо
укладывается в одну траншею и получает UB, либо делает лишнюю докопку и теряет половину. Оценивать
политику раскладки на двух открытых тестах нельзя — нужен набор похожих задач с разными seed.

  python scripts/gen_tiny.py --count 40 --seed 5 --out build/tiny
  python scripts/gen_tiny.py --count 40 --cap 20 --out build/tiny20
"""
import argparse
import os
import random


def gen(rng, cap, maxlen):
    """Случайная последовательность укладок/удалений с пиком ≤ cap (траншея = cap)."""
    cmds = []
    alive = {}          # pid -> length
    live = 0
    pid = 0
    target_lays = rng.randint(4, 14)
    lays = 0
    while lays < target_lays or (alive and rng.random() < 0.5):
        can_lay = live < cap
        if can_lay and (not alive or rng.random() > 0.33):
            space = cap - live
            L = rng.randint(1, max(1, min(maxlen, space)))
            pid += 1
            cmds.append(L)
            alive[pid] = L
            live += L
            lays += 1
        elif alive:
            p = rng.choice(list(alive))
            cmds.append(-p)
            live -= alive.pop(p)
    if not cmds:
        cmds = [1]
    return cmds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--count", type=int, default=40)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--cap", type=int, default=0, help="0 = по кругу 20/40/60")
    ap.add_argument("--maxlen", type=int, default=8)
    ap.add_argument("--out", default="build/tiny")
    args = ap.parse_args()

    caps = [args.cap] if args.cap else [20, 40, 60]
    os.makedirs(args.out, exist_ok=True)
    for i in range(1, args.count + 1):
        rng = random.Random(f"{args.seed}:tiny:{i}")
        cap = caps[(i - 1) % len(caps)]
        cmds = gen(rng, cap, args.maxlen)
        with open(os.path.join(args.out, f"{i}.in"), "w") as f:
            f.write(f"1 {len(cmds)}\n")
            f.write("\n".join(map(str, cmds)))
            f.write("\n0\n")
    print(f"{args.count} tiny tests -> {args.out}")


if __name__ == "__main__":
    main()
