#!/usr/bin/env python3
"""Много маленьких случайных тестов разных классов (для массового fuzz-прогона с --check): build/fuzz2/<seed>/<class>_<i>.in"""
import os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_tests as g1
import gen_tests2 as g2
from gen_tests import write_test
root = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "build", "fuzz2")
seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 5
classes = {}
for name, (fn, sizes) in g1.CLASSES.items(): classes[name] = fn
for name, (fn, sizes) in g2.CLASSES.items(): classes[name] = fn
for seed in range(seeds):
    d = os.path.join(root, str(seed)); os.makedirs(d, exist_ok=True)
    rng0 = random.Random(1000 + seed)
    for name, fn in classes.items():
        for i in range(4):
            n = rng0.choice([40, 150, 600, 2500, 9000, 25000])
            rng = random.Random(f"{seed}:{name}:{i}")
            try:
                cmds = fn(rng, n)
            except Exception as e:
                continue
            if not cmds: continue
            write_test(os.path.join(d, f"{name}_{i}.in"), cmds, 1)
print("done", root)
