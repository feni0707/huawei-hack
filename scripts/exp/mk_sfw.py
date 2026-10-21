#!/usr/bin/env python3
"""mk_sfw.py : solutions/scorefit.cpp -> scripts/exp/sf_w.cpp — веса LW из env (SFW), второй набор LW2 (SFW2) и смесь
по доле «стековых» удалений: a = clamp((lifo8 - BLA) / (BLB - BLA)), веса = a·LW + (1 - a)·LW2. Только для стенда simc."""
src = open('solutions/scorefit.cpp', encoding='utf-8').read()
i = src.index('static const double LW[27] =')
j = src.index('\n', i)
head, body = src[:i], src[j + 1:]
lwline = src[i:j + 1].replace('static const double LW[27] =', 'static double LW[27] =')
body = body.replace('LW[', 'LWb[')
body = body.replace('    static double nbW(int r) {', '    double nbW(int r) const {', 1)
extra = r'''static double LW2[27]; static bool LW2set = false; static double BLA = 0.5, BLB = 0.7;
static void parseVec(const char* e, double* w) { char* q = (char*)e; for (int i = 0; i < 27 && *q; i++) { char* r; double v = strtod(q, &r); if (r == q) break; w[i] = v; q = r; while (*q == ',' || *q == ' ') q++; } }
static void initWeights() {
    static bool done = false; if (done) return; done = true;
    if (const char* e = getenv("SFW")) parseVec(e, LW);
    for (int i = 0; i < 27; i++) LW2[i] = LW[i];
    if (const char* e = getenv("SFW2")) { LW2set = true; parseVec(e, LW2); }
    if (getenv("BLA")) BLA = atof(getenv("BLA"));
    if (getenv("BLB")) BLB = atof(getenv("BLB"));
}
'''
body = body.replace('    ll T = 0;                 // длина траншеи', '    double LWb[27];           // веса для текущей укладки\n    ll T = 0;                 // длина траншеи', 1)
old = '    Alloc() {\n        if (logN[1] == 0) initLogs();'
assert old in body
body = body.replace(old, '    Alloc() {\n        initWeights(); for (int i = 0; i < 27; i++) LWb[i] = LW[i];\n        if (logN[1] == 0) initLogs();', 1)
old = '    void place(int L, ll& outN, ll& outD) {\n        shadowCheck();'
assert old in body
body = body.replace(old, r'''    void place(int L, ll& outN, ll& outD) {
        if (LW2set) {   // смесь: стековые удаления — LW, иначе — LW2
            double a = 1;
            if (lifoTot >= 20) { double h = lifoHit / lifoTot; a = (h - BLA) / (BLB - BLA); if (a < 0) a = 0; if (a > 1) a = 1; }
            for (int i = 0; i < 27; i++) LWb[i] = a * LW[i] + (1 - a) * LW2[i];
        }
        shadowCheck();''', 1)
body = body.replace("        return c + LWb[17] * seal;", "        if (LWb[20] != 0) c += LWb[20] * 100.0 * (double)gs / (double)(T > 0 ? T : 1);   // адрес: <0 — левее (first fit)\n        return c + LWb[17] * seal;", 1)
out = head.replace('#include <cstdio>\n', '#include <cstdio>\n#include <cstdlib>\n', 1) + lwline + extra + body
open(__import__('sys').argv[1] if len(__import__('sys').argv) > 1 else 'scripts/exp/sf_w.cpp', 'w', encoding='utf-8').write(out)
print('ok', out.count('LWb['))
