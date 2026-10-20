// Стенд для быстрых экспериментов: решение подключается через #include (main переименован),
// команды теста подаются напрямую, ответы проверяются ровно как в interactor.cpp, KPI считается так же.
//
//   cl /nologo /O2 /EHsc /std:c++17 /DSIM /DUNIVSRC='"...\solutions\universal.cpp"' /Fe:build\sim2.exe scripts\sim2.cpp
//   build\sim2.exe materials/tests/*.in            — по строке на тест
//   build\sim2.exe --agg build/gen/*              — средний KPI по каждой папке (аргументы — папки)
//   env K="3=1.5 7=0"                             — ручки (массив K[] в решении)
//
// Диагностика: env DIAG=1 печатает число докопок и финальные дыры по размерам; env CUTS=n — KPI, усреднённый
// по n обрезаниям теста в последних 10 % (устойчивость к «где кончился тест»).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include <vector>

typedef long long ll;
double K[64];

std::vector<int> ORD;   // индекс команды удаления по pid (или n) — оракул для диагностики

#define main univ_main_unused
#include UNIVSRC
#undef main

#ifdef USE_ALLOC
typedef Alloc Engine;      // эталонные решения (scorefit.cpp) определяют struct Alloc
#endif
static std::vector<int> readTest(const char* fn, double& w) {
    FILE* f = fopen(fn, "r");
    if (!f) { perror(fn); exit(1); }
    int n;
    if (fscanf(f, "%lf %d", &w, &n) != 2) exit(1);
    std::vector<int> c(n);
    for (int i = 0; i < n; i++) if (fscanf(f, "%d", &c[i]) != 1) exit(1);
    fclose(f);
    return c;
}

struct Result { double kpi = 0; ll T = 0, digs = 0; bool ok = true; const char* why = ""; };

static Result runTest(const std::vector<int>& cmd) {
    Result R;
    Engine* E = new Engine();
    std::vector<int> cells;                 // 0 = свободно, иначе pid
    std::vector<int> plen(1, 0);
    ll occ = 0;
    int n = (int)cmd.size();
    for (int ci = 0; ci < n && R.ok; ci++) {
        int c = cmd[ci];
        if (c == 0) break;
        if (c > 0) {
            ll N = -1, D = -1;
            E->place(c, N, D);
            ll Tnew = (ll)cells.size() + 20 * N;
            if (N < 0) { R.ok = false; R.why = "N<0"; break; }
            if (Tnew != E->T) { R.ok = false; R.why = "T mismatch"; break; }
            if (D < 0 || D + c > Tnew || Tnew > 100000000) { R.ok = false; R.why = "out of trench"; break; }
            cells.resize(Tnew, 0);
            R.digs += N;
            int pid = (int)plen.size();
            plen.push_back(c);
            for (int i = 0; i < c; i++)
                if (cells[D + i]) { R.ok = false; R.why = "collision"; break; }
            if (!R.ok) break;
            for (int i = 0; i < c; i++) cells[D + i] = pid;
            occ += c;
        } else {
            int pid = -c;
            ll d = E->remove(pid);
            if (d < 0 || d >= (ll)cells.size() || cells[d] != pid) { R.ok = false; R.why = "bad position"; break; }
            for (int i = 0; i < plen[pid]; i++) cells[d + i] = 0;
            occ -= plen[pid];
        }
    }
    R.T = E->T;
    R.kpi = (R.ok && R.T > 0) ? (double)occ / R.T : 0.0;
    if (getenv("DIAG") && R.ok) {
        std::vector<int> hist(64, 0);
        for (size_t i = 0; i < cells.size(); i++)
            if (!cells[i]) hist[std::min<size_t>(63, cells.size() - i)]++;   // грубо: не используется
        (void)hist;
    }
    delete E;
    return R;
}

int main(int argc, char** argv) {
#ifndef USE_ALLOC
    for (int i = 0; i < 64; i++) K[i] = KD[i];   // значения по умолчанию из решения
#endif
    if (const char* s = getenv("K")) {
        const char* p = s;
        while (*p) {
            char* q;
            long i = strtol(p, &q, 10);
            if (q == p || *q != '=') break;
            p = q + 1;
            double v = strtod(p, &q);
            if (q == p) break;
            if (i >= 0 && i < 64) K[i] = v;
            p = q;
        }
    }
    bool agg = false;
    (void)agg;
    for (int i = 1; i < argc; i++) {
        double w;
        std::vector<int> cmd = readTest(argv[i], w);
        clock_t t0 = clock();
        Result R = runTest(cmd);
        double cpu = double(clock() - t0) / CLOCKS_PER_SEC;
        printf("%s %g %.6f %lld %lld %.4f %d\n", argv[i], w, R.kpi, R.T, R.digs, cpu, R.ok ? 1 : 0);
        if (!R.ok) fprintf(stderr, "WA %s: %s\n", argv[i], R.why);
    }
    return 0;
}
