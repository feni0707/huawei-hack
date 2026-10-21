#!/usr/bin/env python3
"""mk_sfcls.py [src] [out] : зоны по длинам с классом трубы (ключ зоны = L + 64·класс). Класс задаёт член zCls перед place()
(0 — обычное поведение, побитово как исходник). В стенде (SFSRC) с K[40] > 0 класс берётся из оракула смерти:
класс 1, если труба умрёт не раньше чем через K[40] команд."""
import sys
src = sys.argv[1] if len(sys.argv) > 1 else 'solutions/scorefit_plus.cpp'
out = sys.argv[2] if len(sys.argv) > 2 else 'scripts/exp/sfp_cls.cpp'
s = open(src, encoding='utf-8').read()
def rep(old, new, cnt=1):
    global s
    assert old in s, old[:80]
    s = s.replace(old, new, cnt)
rep('    std::vector<int> slots[MAXL + 1];    // свободные слоты участков длины L (стек)',
    '    std::vector<int> slots[16 * 64];     // свободные слоты участков: ключ = длина + 64·класс (стек), класс < 16\n    int zCls = 0, zKey = 0;              // класс текущей трубы (задаётся снаружи) и ключ её зоны')
rep('    bool zoneOn(int z, bool en) const { return en && z > 0 && cnt[z] > ZSHARE * cntTot; }',
    '    bool zoneOn(int z, bool en) const { return en && z > 0 && cnt[z & 63] > ZSHARE * cntTot; }')
rep('        for (int x = d; x < d + n; x++) { int z = zl.get(x); if (z && z != L && zoneOn(z, en)) return false; }',
    '        for (int x = d; x < d + n; x++) { int z = zl.get(x); if (z && z != zKey && zoneOn(z, en)) return false; }')
rep('            bool ok = mode == 0 ? (zl.get(g) == L || !zoneOn(zl.get(g), en))',
    '            bool ok = mode == 0 ? (zl.get(g) == zKey || !zoneOn(zl.get(g), en))')
rep('        if (zoneOn(L, en) && aliveCntL[L] >= ZMINALIVE) {', '        if (zoneOn(zKey, en) && aliveCntL[L] >= ZMINALIVE) {')
rep('            while (!slots[L].empty()) {\n                int d = slots[L].back(); slots[L].pop_back();\n                if (zl.get(d) != L || occupied(d) || occupied(d + L - 1)) continue;',
    '            while (!slots[zKey].empty()) {\n                int d = slots[zKey].back(); slots[zKey].pop_back();\n                if (zl.get(d) != zKey || occupied(d) || occupied(d + L - 1)) continue;')
rep('                    for (ll x = hg; x < hg + C2; x++) zl.set(x, L);\n                    for (ll x = hg + C2 - L; x > hg; x -= L) slots[L].push_back((int)x);',
    '                    for (ll x = hg; x < hg + C2; x++) zl.set(x, zKey);\n                    for (ll x = hg + C2 - L; x > hg; x -= L) slots[zKey].push_back((int)x);')
rep('            for (ll x = stp; x < stp + C; x++) zl.set(x, L);\n            for (ll x = stp + C - L; x > stp; x -= L) slots[L].push_back((int)x);',
    '            for (ll x = stp; x < stp + C; x++) zl.set(x, zKey);\n            for (ll x = stp + C - L; x > stp; x -= L) slots[zKey].push_back((int)x);')
rep('        if (zl.get(d) == l) slots[l].push_back(d);', '        { int z = zl.get(d); if (z && (z & 63) == l) slots[z].push_back(d); }')
# ключ зоны вычисляется в начале zoneTry
rep('    bool zoneTry(int L, int st, int m, ll& N, int& D) {', '    bool zoneTry(int L, int st, int m, ll& N, int& D) {\n        zKey = L + 64 * zCls;')
# стенд: класс из оракула
orc = r'''
#ifdef SFSRC
extern std::vector<int> ORD; extern long long CUR; extern double K[64];
extern std::vector<unsigned char> gClsTab;   // классы из таблицы (стенд); пусто — оракул K[40]/K[41]
#endif
'''
i = s.index('typedef long long ll;')
s = s[:i] + orc + s[i:]
rep('    void place(int L, ll& outN, ll& outD) {\n', '''    void place(int L, ll& outN, ll& outD) {
#ifdef SFSRC
        if (!gClsTab.empty() && !noStreams) { int pidc = (int)pos.size(); zCls = pidc < (int)gClsTab.size() ? gClsTab[pidc] : 0; }
        else if (K[40] > 0 && !noStreams) { int pidc = (int)pos.size(); long long dt = ORD[pidc] - CUR; zCls = (ORD[pidc] < (1 << 30) && dt >= (long long)K[40]) ? 1 : 0; if (zCls && K[41] > 0) { zCls = 1 + (int)((double)ORD[pidc] / K[41]); if (zCls > 7) zCls = 7; } }
#endif
''')
s = s.replace('shadow = new Alloc(*this); shadow->forceD = -1;', 'shadow = new Alloc(*this); shadow->forceD = -1; shadow->zCls = 0;', 1)
open(out, 'w', encoding='utf-8').write(s)
print('written', out)
