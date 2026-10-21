#!/usr/bin/env python3
"""mk_sfov.py [out] : scorefit_plus + переопределения решений (pid -> номер варианта) + драйвер поиска (-DOVSEARCH).
Сначала применяется mk_sfcf.py (запись вариантов и принудительный выбор), затем добавляется таблица переопределений.
Драйвер: OVSEARCH, режимы (env MODE):
  list : тест + OVFILE -> строки "D pid nopt def" для решений с >= 2 вариантами и итог "B T emax"
  eval : тест + OVFILE + CANDS (файл "pid idx") -> строки "E pid idx T emax" (каждый кандидат добавляется к OVFILE)."""
import subprocess, sys
out = sys.argv[1] if len(sys.argv) > 1 else 'scripts/exp/sfp_ov.cpp'
src = sys.argv[2] if len(sys.argv) > 2 else 'solutions/scorefit_plus.cpp'
subprocess.run([sys.executable, 'scripts/exp/mk_sfcf.py', src, out], check=True)
s = open(out, encoding='utf-8').read()
tbl = r'''
// переопределения решений: отсортированы по pid; номер варианта — индекс в списке вариантов (gOpts) этой укладки
static std::vector<std::pair<int, int>> gOv;
// подпись переопределения (pid, число вариантов, размер выбранного промежутка; хвост = 0): применяется, только если совпала
static std::vector<std::array<int, 3>> gOvSig;
static bool ovSigOk(int pid, const std::vector<CfOpt>& opts, int idx) {
    if (gOvSig.empty()) return true;
    size_t lo = 0, hi = gOvSig.size();
    while (lo < hi) { size_t m = (lo + hi) / 2; if (gOvSig[m][0] < pid) lo = m + 1; else hi = m; }
    if (lo >= gOvSig.size() || gOvSig[lo][0] != pid) return false;
    return (int)opts.size() == gOvSig[lo][1] && (opts[idx].g < 0 ? 0 : opts[idx].s) == gOvSig[lo][2];
}
static int ovFind(int pid) {
    size_t lo = 0, hi = gOv.size();
    while (lo < hi) { size_t m = (lo + hi) / 2; if (gOv[m].first < pid) lo = m + 1; else hi = m; }
    return lo < gOv.size() && gOv[lo].first == pid ? gOv[lo].second : -1;
}
'''
i = s.index('static void cfInit() {')
s = s[:i] + tbl + s[i:]
old = '    void place(int L, ll& outN, ll& outD) {\n'
assert old in s
s = s.replace(old, old + '        if (!gOv.empty() && !noStreams) { int fo = ovFind((int)pos.size()); if (fo >= 0) { gForcePid = (int)pos.size(); gForceIdx = fo; } }\n', 1)
s = s.replace('if (gForcePid == pid && gForceIdx >= 0 && gForceIdx < (int)gOpts.size()) {', 'if (gForcePid == pid && gForceIdx >= 0 && gForceIdx < (int)gOpts.size() && !noStreams && ovSigOk(pid, gOpts, gForceIdx)) {', 1)
s = s.replace('#include <cstdio>\n', '#include <cstdio>\n#include <array>\n', 1)
s = s.replace('static int gForcePid = -1, gForceIdx = -1, gChosen = -1;', 'static int gForcePid = -1, gForceIdx = -1, gChosen = -1;\nstatic bool gRecordForce = false;', 1)
drv = r'''
#ifdef OVSEARCH
#include <cstring>
static std::vector<int> readCmdsO(const char* fn) { FILE* f = fopen(fn, "r"); double w; int n; std::vector<int> c; if (!f || fscanf(f, "%lf %d", &w, &n) != 2) return c; c.resize(n); for (int i = 0; i < n; i++) if (fscanf(f, "%d", &c[i]) != 1) { c.resize(i); break; } fclose(f); return c; }
static std::vector<std::pair<int, int>> readPairs(const char* fn) { std::vector<std::pair<int, int>> v; FILE* f = fn ? fopen(fn, "r") : 0; if (!f) return v; int a, b; while (fscanf(f, "%d %d", &a, &b) == 2) v.push_back({a, b}); fclose(f); return v; }
static void runOv(const std::vector<int>& c, ll& T, ll& emax, bool list) {
    Alloc* A = new Alloc(); emax = 0;
    for (size_t i = 0; i < c.size(); i++) {
        int x = c[i]; if (x == 0) break;
#ifdef OVCLS
        CUR = (long long)i;
#endif
        if (x > 0) {
            int pid = (int)A->pos.size();
            if (list) gRecord = true;
            ll N, D; A->place(x, N, D);
            if (list) { gRecord = false; if (gOpts.size() >= 2) { printf("D %d %zu %d", pid, gOpts.size(), gChosen); for (auto& o : gOpts) printf(" %d", o.g < 0 ? 0 : o.s); printf("\n"); } }
        } else A->remove(-x);
        ll e = A->tailStart(); if (e > emax) emax = e;
    }
    T = A->T; delete A;
}
#ifdef OVCLS
std::vector<int> ORD; long long CUR = 0; double K[64];
std::vector<unsigned char> gClsTab;
#endif
int main(int argc, char** argv) {
    const char* mode = getenv("MODE") ? getenv("MODE") : "list";
    std::vector<int> c = readCmdsO(argv[1]);
#ifdef OVCLS
    ORD.assign(1, 1 << 30); for (size_t q = 0; q < c.size(); q++) { if (c[q] > 0) ORD.push_back(1 << 30); else if (c[q] < 0) ORD[-c[q]] = (int)q; }
    if (getenv("CLSK")) K[40] = atof(getenv("CLSK"));
    if (getenv("CLSB")) K[41] = atof(getenv("CLSB"));
    if (getenv("CLSTAB")) { FILE* f = fopen(getenv("CLSTAB"), "r"); int a, b; while (f && fscanf(f, "%d %d", &a, &b) == 2) { if ((int)gClsTab.size() <= a) gClsTab.resize(a + 1, 0); gClsTab[a] = (unsigned char)b; } if (f) fclose(f); if (gClsTab.empty()) gClsTab.push_back(0); }
    if (getenv("CLSFILE")) {   // классы из файла (по строке «pid» на каждую трубу класса 1): заменяют оракул — ORD подгоняется
        std::vector<char> one(ORD.size(), 0); FILE* f = fopen(getenv("CLSFILE"), "r"); int q; while (f && fscanf(f, "%d", &q) == 1) if (q > 0 && q < (int)one.size()) one[q] = 1; if (f) fclose(f);
        K[40] = 1; std::vector<int> b(1, 0); for (size_t z = 0; z < c.size(); z++) if (c[z] > 0) b.push_back((int)z);
        for (size_t q2 = 1; q2 < ORD.size(); q2++) ORD[q2] = one[q2] ? b[q2] + 1 : (1 << 30);
    }
#endif
    std::vector<std::pair<int, int>> base = readPairs(getenv("OVFILE"));
    std::sort(base.begin(), base.end());
    if (!strcmp(mode, "list")) {
        gOv = base; ll T, em; runOv(c, T, em, true); printf("B %lld %lld\n", T, em);
    } else {
        std::vector<std::pair<int, int>> cand = readPairs(getenv("CANDS"));
        for (auto& cd : cand) {
            gOv = base; bool rep = false;
            for (auto& o : gOv) if (o.first == cd.first) { o.second = cd.second; rep = true; }
            if (!rep) { gOv.push_back(cd); std::sort(gOv.begin(), gOv.end()); }
            gForcePid = -1; gForceIdx = -1;
            ll T, em; runOv(c, T, em, false);
            printf("E %d %d %lld %lld\n", cd.first, cd.second, T, em); fflush(stdout);
        }
    }
    return 0;
}
#define main sf_main_unused_ov
#endif
'''
i = s.index('#ifdef CFDUMP')
s = s[:i] + drv + s[i:]
open(out, 'w', encoding='utf-8').write(s)
print('written', out)
