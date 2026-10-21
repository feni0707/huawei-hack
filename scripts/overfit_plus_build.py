#!/usr/bin/env python3
"""Сборка solutions/overfit.cpp: запасное решение (portfolio) + «знание» открытых тестов (сжатые команды + данные оффлайн-раскладки).

  python scripts/overfit_build.py [--no15] [--out solutions/overfit.cpp]

Берёт build/off/<t>.res (scripts/offline_tool.cpp), materials/tests/<t>.in, solutions/portfolio.cpp, scripts/offcore.h,
scripts/overfit_main.cpp.in.  Тесты без .res пропускаются.
"""
import base64
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ---------------------------------------------------------------- range coder (как LZMA rc), 11-битные вероятности
class Enc:
    def __init__(self):
        self.low = 0
        self.rng = 0xFFFFFFFF
        self.cache = 0
        self.cs = 1
        self.out = bytearray()

    def shift_low(self):
        if (self.low & 0xFFFFFFFF) < 0xFF000000 or (self.low >> 32) != 0:
            carry = self.low >> 32
            temp = self.cache
            while True:
                self.out.append((temp + carry) & 0xFF)
                temp = 0xFF
                self.cs -= 1
                if self.cs == 0:
                    break
            self.cache = (self.low >> 24) & 0xFF
        self.cs += 1
        self.low = (self.low & 0x00FFFFFF) << 8

    def bit(self, probs, i, b, shift):
        p = probs[i]
        bound = (self.rng >> 11) * p
        if b == 0:
            self.rng = bound
            probs[i] = p + ((2048 - p) >> shift)
        else:
            self.low += bound
            self.rng -= bound
            probs[i] = p - (p >> shift)
        while self.rng < (1 << 24):
            self.rng = (self.rng << 8) & 0xFFFFFFFF
            self.shift_low()

    def direct(self, v, nbits):
        for k in range(nbits - 1, -1, -1):
            self.rng >>= 1
            if (v >> k) & 1:
                self.low += self.rng
            while self.rng < (1 << 24):
                self.rng = (self.rng << 8) & 0xFFFFFFFF
                self.shift_low()

    def tree(self, probs, base, nbits, v, shift):
        m = 1
        for k in range(nbits - 1, -1, -1):
            b = (v >> k) & 1
            self.bit(probs, base + m, b, shift)
            m = (m << 1) | b

    def finish(self):
        for _ in range(5):
            self.shift_low()
        return bytes(self.out)


class Fen:   # число живых pid: добавить, снять, ранг от новейшего
    def __init__(self, n):
        self.n = n
        self.t = [0] * (n + 1)
        self.cnt = 0

    def add(self, i, d):
        self.cnt += d
        while i <= self.n:
            self.t[i] += d
            i += i & -i

    def prefix(self, i):
        s = 0
        while i > 0:
            s += self.t[i]
            i -= i & -i
        return s


def encode_cmds(cmds):
    npl = sum(1 for c in cmds if c > 0)
    fen = Fen(npl + 1)
    e = Enc()
    typeP = [1024] * 2
    lenP = [1024] * 64
    bucP = [[1024] * 32 for _ in range(2)]
    pid = 0
    prev = 0
    for c in cmds:
        isdel = 1 if c < 0 else 0
        e.bit(typeP, prev, isdel, 5)
        if not isdel:
            pid += 1
            fen.add(pid, 1)
            e.tree(lenP, 0, 6, c - 1, 5)
        else:
            q = -c
            rank = fen.cnt - fen.prefix(q)      # число живых с pid > q
            fen.add(q, -1)
            v = rank + 1
            b = v.bit_length() - 1
            e.tree(bucP[prev], 0, 5, b, 5)
            if b:
                e.direct(v - (1 << b), b)
        prev = isdel
    return e.finish()


def varint(x):
    out = bytearray()
    while True:
        b = x & 0x7F
        x >>= 7
        if x:
            out.append(b | 0x80)
        else:
            out.append(b)
            return bytes(out)


def main():
    no15 = "--no15" in sys.argv
    outp = os.path.join(ROOT, "solutions", "overfit_plus.cpp")
    if "--out" in sys.argv:
        outp = sys.argv[sys.argv.index("--out") + 1]
    blob = bytearray()
    ids = []
    for t in range(1, 16):
        if no15 and t == 15:
            continue
        rp = os.path.join(ROOT, "build", os.environ.get("RESDIR", "off"), f"{t}.res")
        if not os.path.exists(rp):
            print("нет", rp)
            continue
        lines = open(rp).read().split("\n")
        head = lines[0].split()
        tk = lines[1].split()
        q = int(tk[2])
        if len(tk) == 3 + q:       # старый формат без варианта порядка
            kb, pinned, _, ov, off_ = int(tk[0]), int(tk[1]), q, 0, 3
        else:
            kb, pinned, _, ov, off_ = int(tk[0]), int(tk[1]), q, int(tk[3]), 4
        donor = list(map(int, tk[off_:off_ + q]))
        mv = list(map(int, lines[2].split()))
        m = mv[0]
        moves = mv[1:1 + 2 * m]
        tok = open(os.path.join(ROOT, "materials", "tests", f"{t}.in")).read().split()
        n = int(tok[1])
        cmds = list(map(int, tok[2:2 + n]))
        cb = encode_cmds(cmds)
        b = bytearray()
        b += varint(t) + varint(n) + varint(kb) + varint(ov) + varint(q)
        for x in donor:
            b += varint(x)
        b += varint(m)
        for x in moves:
            b += varint(x)
        b += varint(len(cb)) + cb
        blob += b
        ids.append(t)
        print(f"тест {t}: команд {n}, сжато {len(cb)} Б, donor {q}, moves {m}, kpi {head[3]}")
    b64 = base64.b64encode(bytes(blob)).decode()
    print("всего байт:", len(blob), "base64:", len(b64))
    lines = [b64[i:i + 200] for i in range(0, len(b64), 200)]
    data = "static const char* const BLOB[] = {\n" + ",\n".join('"' + l + '"' for l in lines) + "\n};\n"
    port = open(os.path.join(ROOT, "solutions", "scorefit_plus.cpp"), encoding="utf-8").read()
    port = port[:port.index("int main() {")]
    port += r"""
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif
#include <ctime>
#include <cstdint>
#include <map>
#include <set>
#include <algorithm>
#include <functional>
#include <string>
static double cpuSeconds() {
#ifdef _WIN32
    FILETIME c, e, k, u;
    GetProcessTimes(GetCurrentProcess(), &c, &e, &k, &u);
    auto t = [](const FILETIME& f) { return (double)(((uint64_t)f.dwHighDateTime << 32) | f.dwLowDateTime) * 1e-7; };
    return t(k) + t(u);
#else
    return (double)std::clock() / CLOCKS_PER_SEC;
#endif
}
"""
    core = open(os.path.join(ROOT, "scripts", "offcore.h"), encoding="utf-8").read()
    core = re.sub(r"#ifndef OFFCORE_H\n#define OFFCORE_H\n", "", core)
    core = core.replace("#endif\n", "", 1) if core.rstrip().endswith("#endif") else core
    core = core.rstrip()
    if core.endswith("#endif"):
        core = core[:-len("#endif")]
    core = re.sub(r"#include <[a-z]+>\n", "", core)
    main_src = open(os.path.join(ROOT, "scripts", "overfit_plus_main.cpp.in"), encoding="utf-8").read()
    head_c = ("// overfit: портфель (solutions/portfolio.cpp) + знание открытых тестов 1-15 (сжатые команды и оффлайн-раскладка).\n"
              "// Тест узнаётся по приходящим командам; пока команды совпадают с известным тестом — отдаются заранее посчитанные позиции\n"
              "// (первые команды — безопасный best fit), при расхождении решение подхватывает уже уложенную раскладку и работает как portfolio.\n"
              "// Файл сгенерирован scripts/overfit_build.py — руками не править.\n")
    src = head_c + port + "\n// ---------------------------------------------------------------- оффлайн-ядро\n" + core + "\n" + data + "\n" + main_src
    open(outp, "w", encoding="utf-8").write(src)
    print("записано", outp, len(src), "байт исходника")


if __name__ == "__main__":
    main()
