#!/usr/bin/env python3
"""Собрать пул из подмножества политик другого пула: mkpool.py <src_pool> <name1,name2,...> <dst_pool>"""
import sys
src, names, dst = sys.argv[1], sys.argv[2].split(","), sys.argv[3]
rows = {l.split()[0]: l for l in open(src) if l.strip()}
open(dst, "w").write("".join(rows[n] if rows[n].endswith("\n") else rows[n] + "\n" for n in names))
print(len(names), "->", dst)
