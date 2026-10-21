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
// переопределения позиций: (pid, D, начало и длина свободного промежутка, куда попадает D, в момент укладки — подпись);
// отсортированы по pid. Хвост — промежуток до конца траншеи (длина T - начало), D за концом траншеи — с докопкой.
static std::vector<std::array<int, 4>> gOvD;
static int ovDFind(int pid) {
    size_t lo = 0, hi = gOvD.size();
    while (lo < hi) { size_t m = (lo + hi) / 2; if (gOvD[m][0] < pid) lo = m + 1; else hi = m; }
    return lo < gOvD.size() && gOvD[lo][0] == pid ? (int)lo : -1;
}
'''
i = s.index('static void cfInit() {')
s = s[:i] + tbl + s[i:]
old = '    void place(int L, ll& outN, ll& outD) {\n'
assert old in s
s = s.replace(old, old + '        if (!gOv.empty() && !noStreams) { int fo = ovFind((int)pos.size()); if (fo >= 0) { gForcePid = (int)pos.size(); gForceIdx = fo; } }\n'
              '        if (!gOvD.empty() && !noStreams) { int k = ovDFind((int)pos.size()); if (k >= 0 && dSigOk(gOvD[k], L)) forceD = gOvD[k][1]; }\n', 1)
# проверка подписи позиции (член Alloc) и сброс forceD после укладки
old = '    int forceD = -1;'
assert old in s
s = s.replace(old, '''    // промежуток, куда попадает позиция d при укладке длины L: (начало, длина); хвост — до T (за T — докопка)
    bool gapOf(int d, int L, int& gs, int& gl) {
        int ts = tailStart();
        if (d >= ts) { gs = ts; gl = (int)T - ts; return true; }
        if (occupied(d)) return false;
        gs = d; while (gs > 0 && !occupied(gs - 1)) gs--;
        gl = tag.get(gs);
        return gl > 0 && d + L <= gs + gl;
    }
    bool dSigOk(const std::array<int, 4>& o, int L) { int gs, gl; return gapOf(o[1], L, gs, gl) && gs == o[2] && gl == o[3]; }
''' + old, 1)
old = '            D = forceD;\n'
assert old in s
s = s.replace(old, '            D = forceD; forceD = -1;\n', 1)
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
drvd = r'''
#ifdef DSEARCH
// Поиск переопределений позиций (pid -> D) на одном тесте: случайная труба (pid >= MINPID), случайная позиция среди краёв
// свободных промежутков подходящей длины (и начала хвоста) в момент её укладки, полный прогон; принимается при строгом
// улучшении (T, emax, площадь) — площадь = сумма превышения началом хвоста порога THR по всем командам. Ход «убрать переопределение» — с вер. 0.2.
// env: OUT (файл "pid D gs gl"), DFILE (продолжить), OVFILE (переопределения вариантов: переводятся в позиции),
// MINPID, SECS, SEED; классы — как в OVSEARCH (CLSK/CLSB/CLSTAB).
#include <cstring>
#include <random>
#include <chrono>
static std::vector<int> readCmdsD(const char* fn) { FILE* f = fopen(fn, "r"); double w; int n; std::vector<int> c; if (!f || fscanf(f, "%lf %d", &w, &n) != 2) return c; c.resize(n); for (int i = 0; i < n; i++) if (fscanf(f, "%d", &c[i]) != 1) { c.resize(i); break; } fclose(f); return c; }
// area2 — сумма начала хвоста по всем командам (последний критерий, AREA2: уклон на плато); te — команда, где достигнут emax
static bool gArea2 = false;
static int gProbeD = -1, gProbeGs = 0, gProbeGl = 0; static bool gProbeFail = false;   // CANDD: позиция пробы задана; не поместилась
struct DRes { ll T = 0, emax = 0; double area = 0, area2 = 0; size_t te = 0; bool operator<(const DRes& o) const { if (T != o.T) return T < o.T; if (emax != o.emax) return emax < o.emax; if (area < o.area - 1e-9) return true; if (area > o.area + 1e-9) return false; return gArea2 && area2 < o.area2 - 1e-9; } };
static std::mt19937_64 drng;
static ll gThr = 0; static bool gTrace = false; static long long gVerb = 0;
static double gKn[NKN]; static long long gKnFrom = -1;   // KNFILE/KNFROM: ручки теста и команда, с которой они действуют
static void loadKn() {
    for (int k = 0; k < NKN; k++) gKn[k] = KN_DEF[k];
    if (!getenv("KNFILE")) return;
    FILE* f = fopen(getenv("KNFILE"), "r"); int a; long long q; while (f && fscanf(f, "%d %lld", &a, &q) == 2) if (a >= 0 && a < NKN) gKn[a] = (double)q / 1e4; if (f) fclose(f);
    gKnFrom = getenv("KNFROM") ? atoll(getenv("KNFROM")) : 0;
}
static std::vector<int> gPlan; static std::vector<std::array<int, 4>> gPlanCap; static FILE* gDump = 0;   // PLAN: позиции по pid (-1 — нет); DUMP: "pid D T"   // площадь считается над порогом THR (по умолчанию 0)
// probe >= 0: перед укладкой трубы probe выбрать позицию и поставить переопределение (в gOvD); cap — для записи
// подписей всех укладок с pid из capSet (перевод вариантов в позиции); applied[k] — переопределение k сработало
static DRes runD(const std::vector<int>& c, int probe, std::vector<char>* applied, const std::vector<int>* capSet, std::vector<std::array<int, 4>>* capOut) {
    knSet(KN_DEF);
    Alloc* A = new Alloc(); DRes r;
    if (applied) applied->assign(gOvD.size(), 0);
    for (size_t i = 0; i < c.size(); i++) {
        int x = c[i]; if (x == 0) break;
        if ((long long)i == gKnFrom) knSet(gKn);
#ifdef OVCLS
        CUR = (long long)i;
#endif
        if (x > 0) {
            int pid = (int)A->pos.size();
            if (pid == probe && gProbeD >= 0) {   // заданная позиция (CANDD): подпись — промежуток в этот момент
                int gs, gl;
                if (A->gapOf(gProbeD, x, gs, gl)) {
                    std::array<int, 4> o = {pid, gProbeD, gs, gl}; gProbeGs = gs; gProbeGl = gl;
                    int k = ovDFind(pid);
                    if (k >= 0) gOvD[k] = o; else { gOvD.push_back(o); std::sort(gOvD.begin(), gOvD.end()); if (applied) applied->insert(applied->begin() + ovDFind(pid), 0); }
                } else gProbeFail = true;
            } else if (pid == probe) {
                std::vector<std::array<int, 3>> cand;
                int T = (int)A->T, xx = 0, ts = A->tailStart();
                while (xx < ts) {
                    if (A->occupied(xx)) { int p = A->owner.get(xx); xx = A->pos[p] + A->len[p]; continue; }
                    int gl = A->tag.get(xx); if (gl <= 0) { xx++; continue; }
                    if (gl >= x) { cand.push_back({xx, xx, gl}); if (gl > x) cand.push_back({xx + gl - x, xx, gl}); }
                    xx += gl;
                }
                cand.push_back({ts, ts, T - ts});
                auto cc = cand[drng() % cand.size()];
                std::array<int, 4> o = {pid, cc[0], cc[1], cc[2]};
                int k = ovDFind(pid);
                if (k >= 0) gOvD[k] = o;
                else { gOvD.push_back(o); std::sort(gOvD.begin(), gOvD.end()); if (applied) applied->insert(applied->begin() + ovDFind(pid), 0); }
            }
            if (pid < (int)gPlan.size() && gPlan[pid] >= 0) {
                int gs, gl;
                if (A->gapOf(gPlan[pid], x, gs, gl)) { A->forceD = gPlan[pid]; gPlanCap.push_back({pid, gPlan[pid], gs, gl}); }
                else printf("PLAN: позиция pid %d занята\n", pid);
            }
            ll T0 = A->T; int ts0 = A->tailStart();
            ll N, D; A->place(x, N, D);
            if (gDump) fprintf(gDump, "%d %lld %lld\n", pid, D, A->T);
            if (gTrace && N > 0) {
                int holes = 0, free = 0, xx = 0, mx = 0;
                while (xx < ts0) {
                    if (A->occupied(xx)) { int p = A->owner.get(xx); xx = A->pos[p] + A->len[p]; continue; }
                    int gl = A->tag.get(xx); if (gl <= 0) { xx++; continue; }
                    if (xx < (int)D || xx >= (int)D + x) { free += gl; if (gl >= x) holes++; if (gl > mx) mx = gl; }
                    xx += gl;
                }
                printf("DIG i %zu pid %d L %d T %lld->%lld D %lld ts0 %d holes>=L %d maxhole %d freeIn %d\n", i, pid, x, T0, A->T, D, ts0, holes, mx, free);
            }
            if (applied) { int k = ovDFind(pid); if (k >= 0 && gOvD[k][1] == D) (*applied)[k] = 1; }
            if (capSet && std::binary_search(capSet->begin(), capSet->end(), pid)) {
                int gs, gl;
                if (D >= ts0) { gs = ts0; gl = (int)T0 - ts0; }
                else { gs = (int)D; while (gs > 0 && !A->occupied(gs - 1)) gs--; int e = (int)D + x; while (e < A->T && !A->occupied(e)) e++; gl = e - gs; }
                capOut->push_back({pid, (int)D, gs, gl});
            }
        } else A->remove(-x);
        ll e = A->tailStart(); if (e > r.emax) { r.emax = e; r.te = i; } if (e > gThr) r.area += (double)(e - gThr); r.area2 += (double)e;
    }
    r.T = A->T; delete A; return r;
}
#ifdef OVCLS
std::vector<int> ORD; long long CUR = 0; double K[64];
std::vector<unsigned char> gClsTab;
#endif
static void saveD(const char* fn) { FILE* f = fopen(fn, "w"); for (auto& o : gOvD) fprintf(f, "%d %d %d %d\n", o[0], o[1], o[2], o[3]); fclose(f); }
int main(int argc, char** argv) {
    std::vector<int> c = readCmdsD(argv[1]);
#ifdef OVCLS
    ORD.assign(1, 1 << 30); for (size_t q = 0; q < c.size(); q++) { if (c[q] > 0) ORD.push_back(1 << 30); else if (c[q] < 0) ORD[-c[q]] = (int)q; }
    if (getenv("CLSK")) K[40] = atof(getenv("CLSK"));
    if (getenv("CLSB")) K[41] = atof(getenv("CLSB"));
    if (getenv("CLSTAB")) { FILE* f = fopen(getenv("CLSTAB"), "r"); int a, b; while (f && fscanf(f, "%d %d", &a, &b) == 2) { if ((int)gClsTab.size() <= a) gClsTab.resize(a + 1, 0); gClsTab[a] = (unsigned char)b; } if (f) fclose(f); if (gClsTab.empty()) gClsTab.push_back(0); }
#endif
    const char* out = getenv("OUT") ? getenv("OUT") : "dsearch.out";
    int minpid = getenv("MINPID") ? atoi(getenv("MINPID")) : 1;
    gThr = getenv("THR") ? atoll(getenv("THR")) : 0; gVerb = getenv("VERB") ? atoll(getenv("VERB")) : 0;
    double secs = getenv("SECS") ? atof(getenv("SECS")) : 600;
    loadKn();
    if (getenv("KNDEF")) { for (int k = 0; k < NKN; k++) printf("%d %.6g\n", k, KN_DEF[k]); return 0; }
    if (getenv("KNEVAL")) {
        if (getenv("DFILE")) { FILE* f = fopen(getenv("DFILE"), "r"); int a, b, g1, g2; while (f && fscanf(f, "%d %d %d %d", &a, &b, &g1, &g2) == 4) gOvD.push_back({a, b, g1, g2}); if (f) fclose(f); std::sort(gOvD.begin(), gOvD.end()); }
        DRes r = runD(c, -1, 0, 0, 0);
        printf("K %lld %lld %.0f %.0f\n", r.T, r.emax, r.area, r.area2);
        return 0;
    }
    drng.seed(getenv("SEED") ? atoll(getenv("SEED")) : 1);
    int P = 0; for (int x : c) if (x > 0) P++;
    if (getenv("OVFILE")) {   // переопределения вариантов -> позиции (с подписями этой укладки)
        FILE* f = fopen(getenv("OVFILE"), "r"); int a, b; while (f && fscanf(f, "%d %d", &a, &b) == 2) if (a >= minpid) gOv.push_back({a, b}); if (f) fclose(f);
        std::sort(gOv.begin(), gOv.end());
        std::vector<int> cs; for (auto& o : gOv) cs.push_back(o.first);
        std::vector<std::array<int, 4>> cap;
        DRes r0 = runD(c, -1, 0, &cs, &cap);
        gOv.clear(); gOvD = cap;
        DRes r1 = runD(c, -1, 0, 0, 0);
        printf("перевод: вариантов -> позиций %zu, T %lld/%lld -> %lld/%lld\n", gOvD.size(), r0.T, r0.emax, r1.T, r1.emax);
    }
    if (getenv("CANDD")) {   // оценка кандидатов "pid D" поверх DFILE: строки "E pid D gs gl T emax area area2" (T = -1 — позиция занята)
        if (getenv("DFILE")) { FILE* f = fopen(getenv("DFILE"), "r"); int a, b, g1, g2; while (f && fscanf(f, "%d %d %d %d", &a, &b, &g1, &g2) == 4) gOvD.push_back({a, b, g1, g2}); if (f) fclose(f); std::sort(gOvD.begin(), gOvD.end()); }
        std::vector<std::array<int, 4>> base0 = gOvD;
        FILE* f = fopen(getenv("CANDD"), "r"); int a, b;
        while (f && fscanf(f, "%d %d", &a, &b) == 2) {
            gOvD = base0; gProbeD = b; gProbeFail = false;
            DRes r = runD(c, a, 0, 0, 0);
            if (gProbeFail) printf("E %d %d 0 0 -1 -1 0 0\n", a, b);
            else printf("E %d %d %d %d %lld %lld %.0f %.0f\n", a, b, gProbeGs, gProbeGl, r.T, r.emax, r.area, r.area2);
            fflush(stdout);
        }
        if (f) fclose(f);
        return 0;
    }
    if (getenv("PLAN")) {   // готовая раскладка "pid D" (pid >= MINPID) -> переопределения позиций с подписями
        FILE* f = fopen(getenv("PLAN"), "r"); int a, b; gPlan.assign(P + 2, -1); while (f && fscanf(f, "%d %d", &a, &b) == 2) if (a >= minpid && a <= P) gPlan[a] = b; if (f) fclose(f);
        DRes r0 = runD(c, -1, 0, 0, 0);
        gPlan.clear(); gOvD = gPlanCap;
        DRes r1 = runD(c, -1, 0, 0, 0);
        printf("план: позиций %zu, T %lld/%lld -> %lld/%lld\n", gOvD.size(), r0.T, r0.emax, r1.T, r1.emax);
    }
    if (getenv("DFILE")) { FILE* f = fopen(getenv("DFILE"), "r"); int a, b, g1, g2; while (f && fscanf(f, "%d %d %d %d", &a, &b, &g1, &g2) == 4) gOvD.push_back({a, b, g1, g2}); if (f) fclose(f); std::sort(gOvD.begin(), gOvD.end()); }
    gArea2 = getenv("AREA2") != 0;
    std::vector<int> birthI(1, 0), deathI(1, 1 << 30);
    for (size_t q = 0; q < c.size(); q++) { if (c[q] > 0) { birthI.push_back((int)q); deathI.push_back(1 << 30); } else if (c[q] < 0) deathI[-c[q]] = (int)q; }
    std::vector<int> pool;
    auto mkPool = [&](size_t te) { pool.clear(); for (int q = minpid; q <= P; q++) if (birthI[q] <= (int)te && deathI[q] > (int)te) pool.push_back(q); };
    double pfr = getenv("POOLFR") ? atof(getenv("POOLFR")) : 0.75;
    std::vector<char> ap;
    DRes base = runD(c, -1, &ap, 0, 0);
    { std::vector<std::array<int, 4>> keep; for (size_t k = 0; k < gOvD.size(); k++) if (ap[k]) keep.push_back(gOvD[k]); gOvD = keep; }
    base = runD(c, -1, 0, 0, 0);
    mkPool(base.te);
    printf("start T %lld emax %lld area %.0f overrides %zu placements %d pool %zu (te %zu)\n", base.T, base.emax, base.area, gOvD.size(), P, pool.size(), base.te); fflush(stdout);
    saveD(out);
    auto t0 = std::chrono::steady_clock::now();
    long long it = 0, acc = 0;
    for (;;) {
        double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (el > secs) break;
        it++;
        std::vector<std::array<int, 4>> saved = gOvD;
        DRes r; std::vector<char> a2;
        if (!gOvD.empty() && drng() % 5 == 0) { gOvD.erase(gOvD.begin() + drng() % gOvD.size()); r = runD(c, -1, &a2, 0, 0); }
        else {
            int p = (!pool.empty() && (double)(drng() % 1000) < 1000 * pfr) ? pool[drng() % pool.size()] : minpid + (int)(drng() % (unsigned)(P + 1 - minpid));
            r = runD(c, p, &a2, 0, 0);
        }
        if (gVerb && it <= gVerb) printf("  eval %lld: T %lld emax %lld area %.0f\n", it, r.T, r.emax, r.area);
        if (r < base) {
            std::vector<std::array<int, 4>> keep; for (size_t k = 0; k < gOvD.size(); k++) if (a2[k]) keep.push_back(gOvD[k]);
            bool newT = r.T < base.T;
            gOvD = keep; base = r; acc++; saveD(out); mkPool(base.te);
            if (acc % 20 == 1 || newT) { printf("it %lld %.0fs: T %lld emax %lld area %.0f overrides %zu\n", it, el, base.T, base.emax, base.area, gOvD.size()); fflush(stdout); }
        } else gOvD = saved;
        if (it % 2000 == 0) { printf("it %lld %.0fs: T %lld emax %lld area %.0f overrides %zu acc %lld\n", it, el, base.T, base.emax, base.area, gOvD.size(), acc); fflush(stdout); }
    }
    gTrace = getenv("TRACE") != 0;
    if (getenv("DUMP")) gDump = fopen(getenv("DUMP"), "w");
    DRes chk = runD(c, -1, 0, 0, 0);
    printf("final T %lld emax %lld area %.0f overrides %zu (check T %lld emax %lld) it %lld acc %lld\n", base.T, base.emax, base.area, gOvD.size(), chk.T, chk.emax, it, acc);
    if (gDump) fclose(gDump);
    return 0;
}
#define main sf_main_unused_d
#endif
'''
i = s.index('#ifdef CFDUMP')
s = s[:i] + drv + drvd + s[i:]
open(out, 'w', encoding='utf-8').write(s)
print('written', out)
