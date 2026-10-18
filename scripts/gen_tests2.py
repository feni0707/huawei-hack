#!/usr/bin/env python3
"""Дополнительные классы тестов, откалиброванные по статистике открытых тестов (см. docs/TESTS_ANALYSIS.md, раздел про фазы).

Что выяснилось про открытые тесты (scripts/phase_stats.py, scripts/life_stats*.py):
  * 11 и 15 — ДВЕ ФАЗЫ по n/2 команд. Фаза 1: только короткие трубы 1..7, удаляется ~20 % команд (жертвы — свежие).
    Фаза 2: только длинные 8..M (M=32 в 11, 40 в 15), удаляется ~50 % команд: 82 % — свежие длинные (тяжёлый хвост по возрасту),
    18 % — «обратная заметка» по старым коротким трубам: указатель идёт от верха стопки вниз, удаляя каждую вторую.
    Из-за этого длинные трубы не влезают в дыры 1..7 → ~5 % дыр в финале (потеря KPI у best fit: 0.947 на 15).
  * Длина труб и их срок жизни коррелируют (короткие живут долго) — этим пользуется предсказатель времени жизни.

Классы:
  twophase  — как 11/15 (параметры варьируются: граница короткие/длинные, M, доли удалений, плотность заметки)
  tiny      — крошечные случайные тесты (6–24 команд), как 2/3/4/6 (проверка «точного попадания» и аккуратной упаковки)
  invert    — twophase с перевёрнутой шкалой длин (живучие — длинные): проверка, что решение учится, а не зашито под «короткие»

  python3 scripts/gen_tests2.py --out build/gen2 --seed 1
"""
import argparse
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_tests import Builder, write_test, MAX_LEN  # noqa: E402


def lognormal_rank(rng, mu, sigma):
    return int(math.exp(rng.gauss(mu, sigma)))


def gen_twophase(rng, n, invert=False):
    S = rng.choice([5, 6, 7, 7, 7, 8, 9])                       # граница коротких
    M = rng.choice([24, 28, 32, 36, 40, 40])                     # максимум длинных
    frac1 = rng.uniform(0.45, 0.55)
    p_del1 = rng.uniform(0.15, 0.25)
    p_del2 = rng.uniform(0.45, 0.52)
    old_share = rng.uniform(0.1, 0.25)                           # доля «заметки» по старым коротким в фазе 2
    density = rng.uniform(0.4, 0.6)
    cont = rng.uniform(0.15, 0.4)                                # продолжение серии удалений
    b = Builder(n)
    n1 = int(n * frac1)
    m = sum(min(k, 8) * (cont ** (k - 1)) * (1 - cont) for k in range(1, 60))     # средний размер серии удалений
    q1 = min(0.95, (p_del1 / (1 - p_del1)) / m)
    q2 = min(0.95, (p_del2 / (1 - p_del2)) / m)
    short_stack = []                                             # PID коротких (для заметки)
    ptr = None

    # Для эффективности держим «живые по PID» как список (PID возрастают), удаление по индексу — O(n), но n ≤ 2e5 и редкое.
    live_list = []                                               # живые PID в порядке возрастания
    live_set = set()

    def do_place(L):
        pid = b.place(L)
        live_list.append(pid)
        live_set.add(pid)
        return pid

    def do_remove(pid):
        b.remove(pid)
        live_set.discard(pid)
        # ленивое удаление из live_list: чистим при выборе жертвы
    def pick_from_top(rank):
        # i-й с конца среди живых; чистим мёртвых лениво
        while live_list and live_list[-1] not in live_set:
            live_list.pop()
        cnt = 0
        i = len(live_list) - 1
        while i >= 0:
            if live_list[i] in live_set:
                if cnt == rank:
                    return live_list[i]
                cnt += 1
            i -= 1
        return None

    def pick_from_top_fast(rank):
        # быстрый путь: для малых рангов идём по списку; для больших — случайный индекс (приближённо), без точного подсчёта
        if rank < 400:
            return pick_from_top(rank)
        # приближение: смещаемся от конца на rank позиций (живых доля < 1, поэтому ранг чуть занижается) — для формы хвоста достаточно
        while live_list and live_list[-1] not in live_set:
            live_list.pop()
        i = len(live_list) - 1 - rank
        if i < 0:
            return None
        for j in range(i, -1, -1):
            if live_list[j] in live_set:
                return live_list[j]
        return None

    def short_len():
        return rng.randint(1, S)

    def long_len():
        return rng.randint(S + 1, M)

    # --- фаза 1
    while b.room() and len(b.cmds) < n1:
        if not b.room(2):
            break
        do_place(short_len())
        if live_set and rng.random() < q1:
            k = 1
            while k < 7 and rng.random() < cont:
                k += 1
            for _ in range(k):
                if not b.room() or not live_set:
                    break
                r = int(rng.expovariate(1 / 12.0)) if rng.random() < 0.85 else int(rng.expovariate(1 / 70.0))
                v = pick_from_top_fast(r)
                if v is None:
                    v = live_list[-1] if live_list and live_list[-1] in live_set else None
                if v is None:
                    break
                do_remove(v)
    short_alive = [p for p in live_list if p in live_set]
    ptr = len(short_alive) - 1                                    # указатель заметки (индекс в short_alive сверху вниз)

    # --- фаза 2
    long_pids = []
    while b.room():
        long_pids.append(do_place(long_len()))
        if live_set and rng.random() < q2:
            k = 1
            while k < 10 and rng.random() < cont:
                k += 1
            for _ in range(k):
                if not b.room() or not live_set:
                    break
                if rng.random() < old_share and ptr is not None and ptr >= 0:
                    # заметка: шаг вниз на геометрическое число живых, удаляем
                    while ptr >= 0 and short_alive[ptr] not in live_set:
                        ptr -= 1
                    if ptr < 0:
                        continue
                    skip = 0
                    while rng.random() > density and skip < 20:
                        skip += 1
                    ptr -= skip
                    while ptr >= 0 and short_alive[ptr] not in live_set:
                        ptr -= 1
                    if ptr >= 0:
                        do_remove(short_alive[ptr])
                        ptr -= 1
                else:
                    r = lognormal_rank(rng, 3.7, 2.2)
                    # ранг считаем среди живых длинных: берём из хвоста long_pids
                    while long_pids and long_pids[-1] not in live_set:
                        long_pids.pop()
                    cand = None
                    if long_pids:
                        idx = len(long_pids) - 1 - min(r, len(long_pids) - 1)
                        for j in range(idx, -1, -1):
                            if long_pids[j] in live_set:
                                cand = long_pids[j]
                                break
                    if cand is None:
                        continue
                    do_remove(cand)
    cmds = b.cmds
    if invert:
        # перевернуть шкалу длин: L -> MAX_LEN+1-L (живучими становятся длинные)
        cmds = [(MAX_LEN + 1 - c) if c > 0 else c for c in cmds]
    return cmds


def gen_invert(rng, n):
    return gen_twophase(rng, n, invert=True)



def gen_sorted2(rng, n):
    """Как открытые 5/9: серии одинаковых длин по возрастанию, после серии удаляются трубы серии с фиксированным шагом/смещением
    (в 5/9: шаг 2, смещение 0, в порядке PID). Параметры фиксируются на весь тест."""
    variant = rng.choice(["full", "pow2", "subset", "evens", "pow3"])
    if variant == "full":
        lens = list(range(1, rng.choice([20, 30, 40]) + 1))
    elif variant == "pow2":
        lens = [1, 2, 4, 8, 16, 32]
    elif variant == "pow3":
        lens = [1, 3, 9, 27]
    elif variant == "evens":
        lens = list(range(2, 41, 2))
    else:
        lens = sorted(rng.sample(range(1, MAX_LEN + 1), rng.randint(5, 15)))
    step = rng.choice([2, 2, 2, 3, 4])
    off = rng.choice([0, 0, 0, 1])
    group = max(2 * step, int(n / (1.3 * len(lens))))
    b = Builder(n)
    for L in lens:
        if not b.room(group + group // step + 2):
            break
        pids = [b.place(L) for _ in range(group)]
        for p in pids[off::step]:
            b.remove(p)
    return b.cmds


def gen_tiny(rng, n):
    # n игнорируется: размер выбирается внутри (6..24 команд)
    size = rng.randint(6, 24)
    b = Builder(size)
    maxlen = rng.choice([5, 8, 12, 15])
    p_del = rng.uniform(0.25, 0.45)
    while b.room():
        if b.alive and b.cmds and rng.random() < p_del:
            b.remove_old(rng) if rng.random() < 0.5 else b.remove_recent(rng, 3)
        else:
            b.place(rng.randint(1, maxlen))
    return b.cmds


CLASSES = {
    "twophase": (gen_twophase, [10000, 50000, 200000]),
    "invert": (gen_invert, [10000, 50000]),
    "sorted2": (gen_sorted2, [3000, 7500, 9000]),
    "tiny": (gen_tiny, [0]),
}


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--class", dest="cls", default="all")
    ap.add_argument("--count", type=int, default=0, help="тестов на класс (по умолчанию: 6; tiny — 100)")
    ap.add_argument("--size", type=int, default=0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default=os.path.join(root, "build", "gen2"))
    ap.add_argument("--clean", action="store_true")
    args = ap.parse_args()
    names = list(CLASSES) if args.cls == "all" else [args.cls]
    for name in names:
        fn, sizes = CLASSES[name]
        d = os.path.join(args.out, name)
        os.makedirs(d, exist_ok=True)
        if args.clean:
            for f in os.listdir(d):
                if f.endswith(".in"):
                    os.remove(os.path.join(d, f))
        count = args.count or (100 if name == "tiny" else 6)
        for i in range(1, count + 1):
            size = args.size or sizes[(i - 1) % len(sizes)]
            rng = random.Random(f"{args.seed}:{name}:{i}:{size}")
            cmds = fn(rng, size)
            write_test(os.path.join(d, f"{i}.in"), cmds, 1)
        print(f"{name:10s} {count} tests -> {os.path.relpath(d)}")


if __name__ == "__main__":
    main()
