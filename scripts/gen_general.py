#!/usr/bin/env python3
"""Генератор «общих» тестов: 10 семейств × 10 тестов, НЕ повторяющих структуру открытых тестов.

Цель — проверять стратегию в общем виде, а не на механике 1–15. Каждый тест — случайная комбинация
независимых осей, параметры которых выбираются по seed:
  * распределение длин (равномерное на случайном отрезке, геометрическое, бимодальное, Zipf по случайному
    подмножеству, 1–3 длины, нормальное, «тяжёлые» 30..40, степени двойки);
  * политика удалений (случайная живая, FIFO, LIFO, по возрасту, по длине, самая длинная/короткая);
  * профиль нагрузки — целевое число живых труб во времени (константа, рост, пила, синус, ступеньки).

Семейства (build/general/<family>/NN.in):
  steady     — стационарная нагрузка, случайные удаления
  fifo       — очередь: удаляется самая старая труба (скользящее окно)
  lifo       — стек: удаляется самая новая, пила по числу живых
  age_mix    — смесь «свежих» и «старых» удалений с разными средними возрастами
  by_length  — волны удаления целых классов длин, затем спрос на другие длины
  drift      — распределение длин меняется по фазам, политика удалений — случайная на фазу
  mass       — циклы: накопление → массовое удаление 50–95 % → дозаполнение
  adversary  — фрагментация: заполнить длиной a, выбить часть, затем требовать длины > a; повторять
  extremes   — крайние случаи: крошечные n, только 1 / только 40, без удалений, «уложил-удалил»,
               удалить почти всё в конце, 200 000 команд
  chaos      — каждый сегмент — своя случайная комбинация (длины, политика, профиль)
  (+ семейства по литературе: lifetime … spectrum, robson — нижняя оценка Робсона по Bender et al. 2026)

Формат как в materials/tests: `<weight> <n>`, n команд, `0`. Вес 1. Удаляются только живые PID,
длины 1..40, n ≤ 200 000, в конце остаётся хотя бы одна живая труба.

  python3 scripts/gen_general.py                       # 100 тестов в build/general/, seed 1
  python3 scripts/gen_general.py --seed 5 --out build/general5
  python3 scripts/bench.py build/bin/scorefit --tests build/general/*
"""
import argparse
import math
import os
import random
import shutil

MAX_CMDS = 200_000
MAX_LEN = 40
FAMILY_COUNT = 10


class Builder:
    """Команды + живые PID: общий список, списки по длинам, FIFO-указатель, LIFO-стек (ленивые)."""

    def __init__(self, limit):
        self.cmds = []
        self.limit = limit
        self.pid = 0
        self.alive = []
        self.idx = {}
        self.length = {}
        self.by_len = [[] for _ in range(MAX_LEN + 1)]
        self.lidx = {}
        self.fifo_ptr = 1
        self.stack = []

    def room(self, k=1):
        return len(self.cmds) + k <= self.limit

    def n_alive(self):
        return len(self.alive)

    def place(self, L):
        L = max(1, min(MAX_LEN, int(L)))
        self.pid += 1
        p = self.pid
        self.cmds.append(L)
        self.idx[p] = len(self.alive)
        self.alive.append(p)
        self.length[p] = L
        self.lidx[p] = len(self.by_len[L])
        self.by_len[L].append(p)
        self.stack.append(p)
        return p

    def remove(self, p):
        i = self.idx.pop(p)
        last = self.alive.pop()
        if last != p:
            self.alive[i] = last
            self.idx[last] = i
        lst = self.by_len[self.length[p]]
        j = self.lidx.pop(p)
        last = lst.pop()
        if last != p:
            lst[j] = last
            self.lidx[last] = j
        self.cmds.append(-p)

    # --- выбор жертвы
    def pick_random(self, rng):
        return rng.choice(self.alive)

    def pick_oldest(self, rng):
        while self.fifo_ptr not in self.idx:
            self.fifo_ptr += 1
        return self.fifo_ptr

    def pick_newest(self, rng):
        while self.stack[-1] not in self.idx:
            self.stack.pop()
        return self.stack[-1]

    def pick_recent(self, rng, mean_age):
        for _ in range(10):
            p = self.pid - int(rng.expovariate(1.0 / max(mean_age, 0.5)))
            if p in self.idx:
                return p
        return self.pick_random(rng)

    def pick_old_biased(self, rng, frac):
        """Случайная из самой старой доли `frac` живых (по PID) — приближённо через FIFO-указатель."""
        lo = self.pick_oldest(rng)
        for _ in range(10):
            p = rng.randint(lo, lo + max(1, int((self.pid - lo) * frac)))
            if p in self.idx:
                return p
        return lo

    def pick_by_len(self, rng, lens):
        cand = [L for L in lens if self.by_len[L]]
        if not cand:
            return self.pick_random(rng)
        return rng.choice(self.by_len[rng.choice(cand)])

    def pick_extreme_len(self, rng, longest):
        rng_l = range(MAX_LEN, 0, -1) if longest else range(1, MAX_LEN + 1)
        for L in rng_l:
            if self.by_len[L]:
                return rng.choice(self.by_len[L])
        return self.pick_random(rng)

    def finish(self):
        if not self.alive:
            if len(self.cmds) >= self.limit:
                self.cmds.pop()   # отменяем последнее удаление — эта труба остаётся живой
            else:
                self.cmds.append(1 + (self.pid % MAX_LEN))
        return self.cmds


# ---------------------------------------------------------------- оси: длины

def make_len_dist(rng, kind=None):
    kind = kind or rng.choice(["uniform", "geom", "bimodal", "zipf", "few", "normal", "heavy", "pow2"])
    if kind == "uniform":
        a = rng.randint(1, 30)
        b = rng.randint(a, MAX_LEN)
        return kind, lambda r: r.randint(a, b)
    if kind == "geom":
        mean = rng.uniform(1.5, 12)
        return kind, lambda r: min(MAX_LEN, 1 + int(r.expovariate(1.0 / mean)))
    if kind == "bimodal":
        a = rng.randint(1, 15)
        b = rng.randint(20, MAX_LEN)
        pa = rng.uniform(0.2, 0.8)
        sa = rng.randint(0, 3)
        sb = rng.randint(0, 3)
        return kind, lambda r: (r.randint(max(1, a - sa), a + sa) if r.random() < pa
                                else r.randint(b - sb, min(MAX_LEN, b + sb)))
    if kind == "zipf":
        k = rng.randint(3, 20)
        lens = rng.sample(range(1, MAX_LEN + 1), k)
        s = rng.uniform(0.6, 1.8)
        w = [1.0 / (i + 1) ** s for i in range(k)]
        return kind, lambda r: r.choices(lens, w)[0]
    if kind == "few":
        lens = rng.sample(range(1, MAX_LEN + 1), rng.randint(1, 3))
        w = [rng.random() + 0.1 for _ in lens]
        return kind, lambda r: r.choices(lens, w)[0]
    if kind == "normal":
        mu = rng.uniform(3, 37)
        sd = rng.uniform(1, 10)
        return kind, lambda r: max(1, min(MAX_LEN, int(round(r.gauss(mu, sd)))))
    if kind == "heavy":
        return kind, lambda r: r.randint(30, MAX_LEN) if r.random() < 0.85 else r.randint(1, 29)
    # pow2
    return kind, lambda r: r.choice([1, 2, 4, 8, 16, 32])


# ---------------------------------------------------------------- оси: удаления

def make_policy(rng, kind=None):
    kind = kind or rng.choice(["random", "fifo", "lifo", "recent", "old", "longest", "shortest", "by_len"])
    if kind == "random":
        return kind, lambda b, r: b.pick_random(r)
    if kind == "fifo":
        return kind, lambda b, r: b.pick_oldest(r)
    if kind == "lifo":
        return kind, lambda b, r: b.pick_newest(r)
    if kind == "recent":
        mean = rng.choice([1, 3, 10, 50, 300])
        return kind, lambda b, r: b.pick_recent(r, mean)
    if kind == "old":
        frac = rng.choice([0.05, 0.2, 0.5])
        return kind, lambda b, r: b.pick_old_biased(r, frac)
    if kind == "longest":
        return kind, lambda b, r: b.pick_extreme_len(r, True)
    if kind == "shortest":
        return kind, lambda b, r: b.pick_extreme_len(r, False)
    lens = rng.sample(range(1, MAX_LEN + 1), rng.randint(1, 8))
    return kind, lambda b, r: b.pick_by_len(r, lens)


# ---------------------------------------------------------------- оси: профиль нагрузки

def make_profile(rng, kind=None):
    """Возвращает T(t), t ∈ [0,1] — целевое число живых труб."""
    kind = kind or rng.choice(["const", "grow", "saw", "sine", "steps", "hump"])
    base = int(math.exp(rng.uniform(math.log(5), math.log(20000))))
    if kind == "const":
        return kind, lambda t: base
    if kind == "grow":
        return kind, lambda t: 1 + base * t
    if kind == "saw":
        per = rng.randint(2, 12)
        return kind, lambda t: 1 + base * ((t * per) % 1.0)
    if kind == "sine":
        per = rng.uniform(1, 10)
        return kind, lambda t: 1 + base * (0.5 + 0.5 * math.sin(2 * math.pi * per * t))
    if kind == "steps":
        lv = [int(math.exp(rng.uniform(math.log(2), math.log(base + 2)))) for _ in range(rng.randint(3, 8))]
        return kind, lambda t: lv[min(len(lv) - 1, int(t * len(lv)))]
    return kind, lambda t: 1 + base * math.sin(math.pi * t) ** 2 + base * 0.05


def run_segment(b, rng, n_seg, dist, policy, profile, stiffness=None):
    """n_seg команд: укладка/удаление стохастически тянет число живых к T(t)."""
    k = stiffness if stiffness is not None else rng.uniform(0.5, 4)
    start = len(b.cmds)
    while b.room() and len(b.cmds) - start < n_seg:
        t = (len(b.cmds) - start) / max(1, n_seg)
        T = max(1.0, profile(t))
        a = b.n_alive()
        p_del = 1 / (1 + math.exp(max(-50.0, min(50.0, -k * (a - T) / max(1.0, math.sqrt(T))))))
        if a > 0 and rng.random() < p_del:
            b.remove(policy(b, rng))
        else:
            b.place(dist(rng))


def log_size(rng, lo=1000, hi=MAX_CMDS):
    return int(math.exp(rng.uniform(math.log(lo), math.log(hi))))


# ---------------------------------------------------------------- семейства

def fam_steady(rng, i):
    n = log_size(rng)
    b = Builder(n)
    _, dist = make_len_dist(rng)
    _, pol = make_policy(rng, "random")
    _, prof = make_profile(rng, "const")
    run_segment(b, rng, n, dist, pol, prof)
    return b.finish()


def fam_fifo(rng, i):
    n = log_size(rng)
    b = Builder(n)
    _, dist = make_len_dist(rng)
    _, pol = make_policy(rng, "fifo")
    _, prof = make_profile(rng, rng.choice(["const", "sine", "steps"]))
    run_segment(b, rng, n, dist, pol, prof)
    return b.finish()


def fam_lifo(rng, i):
    n = log_size(rng)
    b = Builder(n)
    _, dist = make_len_dist(rng)
    _, pol = make_policy(rng, "lifo")
    _, prof = make_profile(rng, rng.choice(["saw", "sine", "hump"]))
    run_segment(b, rng, n, dist, pol, prof)
    return b.finish()


def fam_age_mix(rng, i):
    n = log_size(rng)
    b = Builder(n)
    _, dist = make_len_dist(rng)
    pols = [make_policy(rng, k)[1] for k in ("recent", "old", "random")]
    w = [rng.random() for _ in pols]
    _, prof = make_profile(rng)
    run_segment(b, rng, n, dist, lambda bb, r: r.choices(pols, w)[0](bb, r), prof)
    return b.finish()


def fam_by_length(rng, i):
    n = log_size(rng, 5000)
    b = Builder(n)
    waves = rng.randint(2, 6)
    for _ in range(waves):
        # набор живых с одним распределением, затем выбиваем выбранные длины и кладём другое распределение
        _, dist1 = make_len_dist(rng)
        fill = n // (waves * 3)
        for _ in range(fill):
            if not b.room():
                break
            b.place(dist1(rng))
        present = [L for L in range(1, MAX_LEN + 1) if b.by_len[L]]
        if not present:
            continue
        victims = rng.sample(present, max(1, len(present) // rng.randint(2, 4)))
        frac = rng.uniform(0.4, 1.0)
        for L in victims:
            for p in list(b.by_len[L]):
                if b.room() and rng.random() < frac:
                    b.remove(p)
        _, dist2 = make_len_dist(rng)
        _, prof = make_profile(rng, "const")
        base = max(5, b.n_alive())
        run_segment(b, rng, n // (waves * 3), dist2, make_policy(rng)[1], lambda t: base)
    return b.finish()


def fam_drift(rng, i):
    n = log_size(rng)
    b = Builder(n)
    phases = rng.randint(3, 10)
    _, prof = make_profile(rng)
    for ph in range(phases):
        _, dist = make_len_dist(rng)
        _, pol = make_policy(rng)
        off, span = ph / phases, 1 / phases
        run_segment(b, rng, n // phases, dist, pol, lambda t, o=off, s=span: prof(o + s * t))
    return b.finish()


def fam_mass(rng, i):
    n = log_size(rng, 3000)
    b = Builder(n)
    _, dist = make_len_dist(rng)
    cycles = rng.randint(2, 8)
    for _ in range(cycles):
        grow = n // (cycles * 2)
        for _ in range(grow):
            if not b.room():
                break
            if b.alive and rng.random() < 0.15:
                b.remove(b.pick_random(rng))
            else:
                b.place(dist(rng))
        _, pol = make_policy(rng)
        kill = int(b.n_alive() * rng.uniform(0.5, 0.95))
        for _ in range(kill):
            if not b.room() or b.n_alive() <= 1:
                break
            b.remove(pol(b, rng))
        if rng.random() < 0.5:
            _, dist = make_len_dist(rng)
    return b.finish()


def fam_adversary(rng, i):
    n = log_size(rng, 3000)
    b = Builder(n)
    mode = rng.choice(["every_other", "random_holes", "grow_len", "prime_gaps"])
    a = rng.randint(1, 8)
    while b.room(10):
        if mode == "grow_len":
            L1, L2 = a, min(MAX_LEN, a + rng.randint(1, 3))
        elif mode == "prime_gaps":
            L1 = rng.choice([2, 3, 5, 7, 11, 13])
            L2 = rng.choice([L for L in (2, 3, 5, 7, 11, 13, 17, 19, 23) if L != L1])
        else:
            L1 = rng.randint(1, 20)
            L2 = rng.randint(L1 + 1, min(MAX_LEN, 2 * L1 + 1))
        cnt = max(4, min(rng.randint(20, 3000), (b.limit - len(b.cmds)) // 3))
        pids = []
        for _ in range(cnt):
            if not b.room():
                break
            pids.append(b.place(L1))
        if mode in ("every_other", "grow_len", "prime_gaps"):
            vict = pids[rng.randint(0, 1)::2]
        else:
            vict = [p for p in pids if rng.random() < rng.uniform(0.3, 0.7)]
        for p in vict:
            if b.room():
                b.remove(p)
        for _ in range(len(vict) // rng.randint(1, 3)):
            if b.room():
                b.place(L2)
        a = a + rng.randint(1, 4) if a < 30 else rng.randint(1, 8)
    return b.finish()


def fam_extremes(rng, i):
    kinds = ["tiny", "all1", "all40", "insert_only", "place_remove", "delete_most",
             "max_n_long", "max_n_mixed", "single_len_churn", "one_survivor"]
    kind = kinds[i % len(kinds)]
    if kind == "tiny":
        n = rng.randint(1, 12)
        b = Builder(n)
        while b.room():
            if b.alive and rng.random() < 0.4:
                b.remove(b.pick_random(rng))
            else:
                b.place(rng.randint(1, MAX_LEN))
        return b.finish()
    if kind in ("all1", "all40"):
        L = 1 if kind == "all1" else MAX_LEN
        n = log_size(rng, 10000)
        b = Builder(n)
        run_segment(b, rng, n, lambda r: L, make_policy(rng)[1], make_profile(rng)[1])
        return b.finish()
    if kind == "insert_only":
        n = log_size(rng, 10000)
        b = Builder(n)
        _, dist = make_len_dist(rng)
        while b.room():
            b.place(dist(rng))
        return b.finish()
    if kind == "place_remove":
        n = log_size(rng, 10000)
        b = Builder(n)
        _, dist = make_len_dist(rng)
        keep = rng.uniform(0.0, 0.1)
        while b.room(2):
            p = b.place(dist(rng))
            if rng.random() >= keep:
                b.remove(p)
        return b.finish()
    if kind == "delete_most":
        n = log_size(rng, 10000)
        b = Builder(n)
        _, dist = make_len_dist(rng)
        while b.room() and len(b.cmds) < n * 0.55:
            b.place(dist(rng))
        _, pol = make_policy(rng)
        while b.room() and b.n_alive() > 1:
            b.remove(pol(b, rng))
        return b.finish()
    if kind == "max_n_long":
        b = Builder(MAX_CMDS)
        run_segment(b, rng, MAX_CMDS, lambda r: r.randint(35, MAX_LEN),
                    make_policy(rng, "random")[1], lambda t: 50000)
        return b.finish()
    if kind == "max_n_mixed":
        b = Builder(MAX_CMDS)
        run_segment(b, rng, MAX_CMDS, lambda r: r.randint(1, MAX_LEN),
                    make_policy(rng, "random")[1], make_profile(rng, "sine")[1])
        return b.finish()
    if kind == "single_len_churn":
        n = log_size(rng, 20000)
        b = Builder(n)
        L = rng.randint(1, MAX_LEN)
        run_segment(b, rng, n, lambda r: L, make_policy(rng, "random")[1],
                    lambda t: 1000, stiffness=0.3)
        return b.finish()
    # one_survivor: всё, кроме случайной ранней трубы, удаляется, потом немного укладок
    n = log_size(rng, 5000)
    b = Builder(n)
    _, dist = make_len_dist(rng)
    while b.room() and len(b.cmds) < n * 0.45:
        b.place(dist(rng))
    survivor = rng.choice(b.alive)
    for p in list(b.alive):
        if p != survivor and b.room():
            b.remove(p)
    while b.room():
        b.place(dist(rng))
    return b.finish()


def fam_chaos(rng, i):
    n = log_size(rng)
    b = Builder(n)
    segs = rng.randint(2, 12)
    for _ in range(segs):
        _, dist = make_len_dist(rng)
        _, pol = make_policy(rng)
        _, prof = make_profile(rng)
        run_segment(b, rng, int(n * rng.uniform(0.3, 1.7) / segs), dist, pol, prof)
    # остаток — последней комбинацией
    run_segment(b, rng, n, make_len_dist(rng)[1], make_policy(rng)[1], make_profile(rng)[1])
    return b.finish()


# ---------------------------------------------------------------- семейства из литературы (docs/RESEARCH_PROMPT_TESTS.md)
# Время в моделях со временем жизни — номер укладки; удаление выдаётся, как только наступил срок смерти трубы.

def make_lifetime(rng, kind=None):
    """Распределение времени жизни (в укладках)."""
    kind = kind or rng.choice(["exp", "pareto", "bimodal", "const", "uniform", "lognormal"])
    mean = math.exp(rng.uniform(math.log(5), math.log(20000)))
    if kind == "exp":                      # Coffman–Leighton / Knuth: пуассоновский поток, M/M/∞
        return kind, lambda r: r.expovariate(1 / mean)
    if kind == "pareto":                   # тяжёлый хвост: большинство живёт мало, единицы — вечно
        a = rng.uniform(1.1, 2.5)
        xm = mean * (a - 1) / a
        return kind, lambda r: xm * r.paretovariate(a)
    if kind == "bimodal":                  # generational hypothesis: «умирают молодыми» + долгожители
        p = rng.uniform(0.6, 0.97)
        short = rng.uniform(1, 30)
        return kind, lambda r: r.expovariate(1 / short) if r.random() < p else r.expovariate(1 / (mean * 5))
    if kind == "const":                    # одинаковое время жизни: смерти строго в порядке рождения
        return kind, lambda r: mean
    if kind == "uniform":
        return kind, lambda r: r.uniform(0, 2 * mean)
    sd = rng.uniform(0.5, 2)
    return kind, lambda r: math.exp(r.gauss(math.log(mean) - sd * sd / 2, sd))


def run_lifetimes(b, rng, n, size_of, life_of, keep_end=0.0):
    """Укладки по одной на тик; смерти — по сроку. size_of(r) -> (L, tag); life_of(r, L, tag) -> время жизни.
    keep_end — доля команд в конце, когда новые трубы получают «вечную» жизнь (их не удалят)."""
    import heapq
    heap = []
    t = 0
    while b.room():
        while heap and heap[0][0] <= t and b.room():
            _, p = heapq.heappop(heap)
            if p in b.idx:
                b.remove(p)
        if not b.room():
            break
        L, tag = size_of(rng)
        p = b.place(L)
        if len(b.cmds) < n * (1 - keep_end):
            heapq.heappush(heap, (t + max(0.5, life_of(rng, L, tag)), p))
        t += 1


def fam_lifetime(rng, i):
    """M/G/∞: независимые времена жизни разных форм; иногда время жизни коррелирует с длиной (Barrett–Zorn)."""
    n = log_size(rng)
    b = Builder(n)
    _, dist = make_len_dist(rng)
    _, life = make_lifetime(rng)
    corr = rng.choice([0, 0, 1, -1])       # 1 — длинные живут дольше, −1 — короче
    f = (lambda r, L, tag: life(r) * (L / 20.5) ** (2 * corr)) if corr else (lambda r, L, tag: life(r))
    run_lifetimes(b, rng, n, lambda r: (dist(r), 0), f)
    return b.finish()


def fam_kv_update(rng, i):
    """Хранилище ключ–значение / журнал (SepBIT, Facebook memcached): обновление ключа удаляет старую версию и
    кладёт новую (длина та же или чуть меняется). Ключи — Zipf: горячие перезаписываются постоянно."""
    n = log_size(rng)
    b = Builder(n)
    nkeys = int(math.exp(rng.uniform(math.log(20), math.log(30000))))
    s = rng.uniform(0.5, 1.4)
    _, dist = make_len_dist(rng)
    jitter = rng.choice([0, 0, 1, 3])      # насколько меняется длина значения при обновлении
    p_new = rng.uniform(0.0, 0.3)          # доля «вставок нового ключа» (рост) против обновлений
    p_del = rng.uniform(0.0, 0.15)         # доля удалений ключа без перезаписи
    cdf, acc = [], 0.0
    for k in range(nkeys):
        acc += 1 / (k + 1) ** s
        cdf.append(acc)
    perm = list(range(nkeys)); rng.shuffle(perm)
    import bisect
    ver, size = {}, {}                      # ключ -> живой PID, длина
    while b.room(2):
        k = perm[bisect.bisect_left(cdf, rng.random() * acc)]
        if k in ver:
            u = rng.random()
            if u < p_del:
                b.remove(ver.pop(k))
                continue
            b.remove(ver[k])
            L = max(1, min(MAX_LEN, size[k] + (rng.randint(-jitter, jitter) if jitter else 0)))
        else:
            if ver and rng.random() > p_new and len(ver) > 5:
                continue
            L = dist(rng)
        size[k] = L
        ver[k] = b.place(L)
    return b.finish()


def fam_cohorts(rng, i):
    """Wilson/Johnstone: объекты, рождённые вместе, умирают вместе («peaks»). Пачки разной природы
    перемежаются, пачка умирает целиком (или почти) через случайное время; часть пачек — «плато» (вечные)."""
    n = log_size(rng, 2000)
    b = Builder(n)
    import heapq
    heap, t = [], 0
    p_plateau = rng.uniform(0.0, 0.4)
    survive = rng.choice([0.0, 0.0, 0.05, 0.2])
    interleave = rng.random() < 0.5         # пачки укладываются вперемешку (2–3 одновременно)
    while b.room():
        while heap and heap[0][0] <= t and b.room():
            _, batch = heapq.heappop(heap)
            order = list(batch)
            if rng.random() < 0.5:
                rng.shuffle(order)
            for p in order:
                if p in b.idx and b.room() and rng.random() >= survive:
                    b.remove(p)
        k = int(math.exp(rng.uniform(math.log(3), math.log(3000))))
        nb = rng.randint(2, 3) if interleave else 1
        batches = [[] for _ in range(nb)]
        dists = [make_len_dist(rng)[1] for _ in range(nb)]
        for _ in range(k):
            if not b.room():
                break
            j = rng.randrange(nb)
            batches[j].append(b.place(dists[j](rng)))
            t += 1
        for bt in batches:
            if bt and rng.random() >= p_plateau:
                heapq.heappush(heap, (t + rng.expovariate(1 / max(1, k * rng.uniform(0.5, 8))), bt))
    return b.finish()


def fam_program(rng, i):
    """Синтетическая «программа»: несколько сайтов выделения (своя длина и своё время жизни у каждого,
    Barrett–Zorn), сайты включаются фазами — ramp (копит), peak (короткоживущие), plateau (долгоживущие)."""
    n = log_size(rng)
    b = Builder(n)
    nsites = rng.randint(3, 12)
    sites = []
    for _ in range(nsites):
        lens = rng.sample(range(1, MAX_LEN + 1), rng.randint(1, 3))
        kind = rng.choice(["short", "short", "medium", "long", "forever"])
        mean = {"short": rng.uniform(1, 20), "medium": rng.uniform(50, 1000),
                "long": rng.uniform(2000, 50000), "forever": 1e12}[kind]
        sites.append((lens, mean))
    nph = rng.randint(1, 8)
    import heapq
    heap, t = [], 0
    for ph in range(nph):
        active = rng.sample(range(nsites), rng.randint(1, min(4, nsites)))
        w = [rng.random() + 0.05 for _ in active]
        end = len(b.cmds) + n // nph
        while b.room() and len(b.cmds) < end:
            while heap and heap[0][0] <= t and b.room():
                _, p = heapq.heappop(heap)
                if p in b.idx:
                    b.remove(p)
            if not b.room():
                break
            s = rng.choices(active, w)[0]
            lens, mean = sites[s]
            p = b.place(rng.choice(lens))
            heapq.heappush(heap, (t + rng.expovariate(1 / mean), p))
            t += 1
    return b.finish()


def fam_generational(rng, i):
    """Generational hypothesis: почти все умирают очень молодыми (часто сразу), выжившие — навсегда;
    доля выживших и их длины меняются (как старение кучи / долгоживущий кэш)."""
    n = log_size(rng)
    b = Builder(n)
    _, d_short = make_len_dist(rng)
    _, d_long = make_len_dist(rng) if rng.random() < 0.6 else (None, d_short)
    p_surv = rng.uniform(0.02, 0.35)
    young = rng.uniform(1, 15)
    import heapq
    heap, t = [], 0
    while b.room():
        while heap and heap[0][0] <= t and b.room():
            _, p = heapq.heappop(heap)
            if p in b.idx:
                b.remove(p)
        if not b.room():
            break
        if rng.random() < p_surv:
            b.place(d_long(rng))
        else:
            p = b.place(d_short(rng))
            heapq.heappush(heap, (t + rng.expovariate(1 / young), p))
        t += 1
    return b.finish()


def fam_slivers(rng, i):
    """Анти-best-fit (по мотивам Robson): спрос на длину h−d как раз тогда, когда освободились дыры h, — каждое
    точное по размеру решение оставляет щель d, щели копятся. Затем запросы длиннее любой дыры."""
    n = log_size(rng, 3000)
    b = Builder(n)
    d = rng.randint(1, 3)
    while b.room(10):
        h = rng.randint(d + 1, MAX_LEN)
        sep = rng.randint(1, MAX_LEN)
        cnt = max(2, min(rng.randint(10, 2000), (b.limit - len(b.cmds)) // 4))
        holes = []
        for _ in range(cnt):
            if not b.room(2):
                break
            holes.append(b.place(h))
            b.place(sep)
        for p in holes:
            if b.room() and rng.random() < 0.9:
                b.remove(p)
        for _ in range(len(holes)):
            if b.room():
                b.place(h - d if rng.random() < 0.8 else rng.randint(1, MAX_LEN))
        # иногда — удалить разделители, чтобы проверить слияние щелей
        if rng.random() < 0.3:
            for p in list(b.alive):
                if b.room() and b.length[p] == sep and rng.random() < 0.5:
                    b.remove(p)
    return b.finish()


def fam_size_shift(rng, i):
    """Нестационарность: распределение длин плавно смещается (короткие → длинные или наоборот), удаления
    случайные/по возрасту; старые дыры чужого размера становятся бесполезны, если не слились."""
    n = log_size(rng)
    b = Builder(n)
    up = rng.random() < 0.75
    lo1, hi1 = (1, rng.randint(3, 10)) if up else (rng.randint(20, 35), MAX_LEN)
    lo2, hi2 = (rng.randint(15, 30), MAX_LEN) if up else (1, rng.randint(3, 10))
    _, pol = make_policy(rng)
    _, prof = make_profile(rng)
    k = rng.uniform(0.5, 4)

    def dist_at(r, t):
        if r.random() < t:
            return r.randint(lo2, hi2)
        return r.randint(lo1, hi1)
    start = len(b.cmds)
    while b.room():
        t = (len(b.cmds) - start) / n
        T = max(1.0, prof(t))
        a = b.n_alive()
        p_del = 1 / (1 + math.exp(max(-50.0, min(50.0, -k * (a - T) / max(1.0, math.sqrt(T))))))
        if a > 0 and rng.random() < p_del:
            b.remove(pol(b, rng))
        else:
            b.place(dist_at(rng, t))
    return b.finish()


def fam_spectrum(rng, i):
    """Выделение спектра (EON) / смежных узлов: запросы размером из малого набора (1, 2, 4, 8 … или 3 типа
    сервисов), пуассоновский приход, экспоненциальное удержание, нагрузка в Эрлангах задаёт число живых;
    в середине — скачок нагрузки (час пик)."""
    n = log_size(rng)
    b = Builder(n)
    sizes = rng.choice([[1, 2, 4, 8, 16, 32], [3, 6, 12, 24], rng.sample(range(1, MAX_LEN + 1), 3),
                        [2, 5, 10, 20, 40]])
    w = [rng.random() + 0.1 for _ in sizes]
    erl = math.exp(rng.uniform(math.log(10), math.log(20000)))
    peak = rng.uniform(1, 4)
    pk0, pk1 = sorted([rng.random(), rng.random()])
    while b.room():
        t = len(b.cmds) / n
        T = erl * (peak if pk0 <= t <= pk1 else 1)
        a = b.n_alive()
        # M/M/∞: доля событий-уходов = a / (a + T)
        if a and rng.random() < a / (a + T):
            b.remove(b.pick_random(rng))
        else:
            b.place(rng.choices(sizes, w)[0])
    return b.finish()


def fam_robson(rng, i):
    """Нижняя оценка Робсона в рандомизированной форме Bender et al. 2026 (arXiv 2608.28462, §3.1): раунды
    r = 0..R, в каждом сначала удаляется случайная доля q живых труб КАЖДОГО размера, затем укладывается объём
    ≈ M размером s_r (M/s_r труб). Размеры растут (1, 2, 4 … 32 или другая цепочка ×2..×3), поэтому старые
    мелкие трубы, умирая вразброс, оставляют дыры, в которые следующие размеры не влезают. Против этого
    теста любой онлайн-аллокатор теряет Ω(log M) — цель не выиграть, а не провалиться сильнее best fit.
    Серий S (весь подъём размеров повторяется); иногда размеры идут вниз или раунд немного «шумит» по длине."""
    n = log_size(rng, 3000)
    b = Builder(n)
    kind = i % 5
    if kind == 0:
        sizes = [1, 2, 4, 8, 16, 32]
    elif kind == 1:
        sizes = [1, 3, 9, 27]
    elif kind == 2:
        sizes = [5, 10, 20, 40]
    elif kind == 3:                                   # случайная цепочка с шагом ×2..×3
        s, sizes = rng.randint(1, 3), []
        while s <= MAX_LEN:
            sizes.append(s)
            s = max(s + 1, int(round(s * rng.uniform(2, 3))))
    else:                                             # плотная цепочка: 1, 2, 3, 5, 8, 13 … (рост ×1.6)
        s, sizes = 1, []
        while s <= MAX_LEN:
            sizes.append(s)
            s = max(s + 1, int(round(s * 1.6)))
    q = rng.uniform(0.05, 0.3)
    series = rng.choice([1, 2, 3, 5])
    down = rng.random() < 0.2
    jitter = rng.random() < 0.3
    # команд на серию ≈ Σ_r M/s_r · (1 + q·раундов) — подбираем M под бюджет
    per_m = sum(1.0 / s for s in sizes) * (1 + q * len(sizes))
    M = max(len(sizes) * MAX_LEN, int(n / series / per_m))
    while b.room(2):
        order = sizes[::-1] if down else sizes
        for s in order:
            for L in range(1, MAX_LEN + 1):
                lst = b.by_len[L]
                if not lst:
                    continue
                for p in rng.sample(lst, int(len(lst) * q + rng.random())):
                    if b.room(2):
                        b.remove(p)
            for _ in range(max(1, M // s)):
                if not b.room(2):
                    break
                L = s
                if jitter and rng.random() < 0.1:
                    L = max(1, min(MAX_LEN, s + rng.choice([-1, 1])))
                b.place(L)
            if not b.room(2):
                break
    return b.finish()


FAMILIES = {
    "steady": fam_steady, "fifo": fam_fifo, "lifo": fam_lifo, "age_mix": fam_age_mix,
    "by_length": fam_by_length, "drift": fam_drift, "mass": fam_mass,
    "adversary": fam_adversary, "extremes": fam_extremes, "chaos": fam_chaos,
    # по литературе
    "lifetime": fam_lifetime, "kv_update": fam_kv_update, "cohorts": fam_cohorts, "program": fam_program,
    "generational": fam_generational, "slivers": fam_slivers, "size_shift": fam_size_shift,
    "spectrum": fam_spectrum, "robson": fam_robson,
}


def validate(cmds):
    assert 1 <= len(cmds) <= MAX_CMDS, len(cmds)
    alive, pid = set(), 0
    for c in cmds:
        if c > 0:
            assert c <= MAX_LEN
            pid += 1
            alive.add(pid)
        else:
            assert -c in alive, c
            alive.remove(-c)
    assert alive, "no pipes left at the end"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--count", type=int, default=FAMILY_COUNT, help="тестов на семейство")
    ap.add_argument("--family", choices=list(FAMILIES), nargs="*", help="по умолчанию все")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "build", "general"))
    ap.add_argument("--clean", action="store_true", help="удалить --out перед генерацией")
    args = ap.parse_args()
    if args.clean and os.path.isdir(args.out):
        shutil.rmtree(args.out)
    for fi, name in enumerate(args.family or FAMILIES):
        d = os.path.join(args.out, name)
        os.makedirs(d, exist_ok=True)
        for i in range(args.count):
            rng = random.Random(f"{args.seed}/{name}/{i}")
            cmds = FAMILIES[name](rng, i)
            validate(cmds)
            with open(os.path.join(d, f"{i + 1:02d}.in"), "w") as f:
                f.write(f"1 {len(cmds)}\n")
                f.write("\n".join(map(str, cmds)))
                f.write("\n0\n")
        print(f"{name:10s} {args.count} tests -> {d}")


if __name__ == "__main__":
    main()
