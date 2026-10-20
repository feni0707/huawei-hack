#!/usr/bin/env python3
"""Граничные и «шумовые» тесты для проверки корректности/времени/памяти: build/fuzz/*.in"""
import os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_tests import write_test, Builder

out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "build", "fuzz")
os.makedirs(out, exist_ok=True)
rng = random.Random(12345)

def save(name, cmds, w=1):
    write_test(os.path.join(out, name + ".in"), cmds, w)

save("single", [7])
save("single_del", [1, -1])
save("max_len_200k", [40] * 200000)
save("len1_200k", [1] * 200000)
# все укладки, потом удалить всё кроме каждой 3-й
b = Builder(200000)
n1 = 100000
pids = [b.place(rng.randint(1, 40)) for _ in range(n1)]
for p in pids:
    if p % 3 and b.room(): b.remove(p)
while b.room(): b.place(rng.randint(1, 40))
save("bulk_delete", b.cmds)
# шум: место/удаление случайной живой, длины равномерно
for name, pdel, lens in (("noise_50", 0.5, (1, 40)), ("noise_30", 0.3, (1, 40)), ("noise_short", 0.45, (1, 4)), ("noise_long", 0.45, (30, 40))):
    b = Builder(200000)
    while b.room():
        if b.alive and rng.random() < pdel: b.remove_old(rng)
        else: b.place(rng.randint(*lens))
    save(name, b.cmds)
# лесенка: места растут, удаления свежих
b = Builder(100000)
while b.room(2):
    L = rng.randint(1, 40)
    b.place(L)
    if rng.random() < 0.9: b.remove(b.pid)
save("place_del_last", b.cmds)
# после полного опустошения (все удалены) — снова укладки
b = Builder(50000)
for _ in range(2000): b.place(rng.randint(1, 40))
for p in list(b.alive): b.remove(p)
while b.room(): b.place(rng.randint(1, 40))
save("empty_then_refill", b.cmds)
print("written to", out)
