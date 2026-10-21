// Оффлайн-стенд для любого решения-класса Alloc (scorefit и копии): решение подключается через #include
// с переименованным main, команды теста подаются напрямую, ответы проверяются как в интеракторе.
//   c++ -std=c++17 -O2 -DSFSRC='"'$PWD/solutions/scorefit.cpp'"' -o build/simc scripts/simc.cpp
//   build/simc materials/tests/*.in          -> "file w=.. kpi=.. T=.." (WA в конце строки при ошибке)
// Ручки экспериментов: в копии решения объявить `extern double K[64];`, задавать env K="i=v i=v".
// Диагностика (env): SFCUTS=20 — KPI усредняется по 20 концам в последних 10 % (меньше шума «где кончился тест»);
// DIAG=1 — свободные клетки по длине бывшей трубы; HOLES=файл — финальные дыры (размер, возраст, длина, занята ли);
// DIGS=1 — каждая докопка со списком дыр; DUMPAT=k DUMPF=файл — раскладка клеток перед командой k;
// NEVSEG=1 (с DIAG) — отрезки ни разу не занятых клеток. ORD/CUR — оракул смерти (время удаления по PID).
// С -DLWMUT веса LW берутся из env SFW (в копии решения убрать const у LW).
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
std::vector<int> ORD; long long CUR = 0, NCMD = 0;
std::vector<unsigned char> gClsTab;   // классы труб для стенда (scripts/exp/mk_sfcls.py)
double K[64]; // ручки экспериментов: env K="i=v i=v"
#define main sf_main_unused
#include SFSRC
#undef main
#include <vector>
static std::vector<int> readTest(const char* fn, double& w) {
    FILE* f = fopen(fn, "r"); if (!f) { perror(fn); exit(1); }
    int n; if (fscanf(f, "%lf %d", &w, &n) != 2) exit(1);
    std::vector<int> c(n);
    for (int i = 0; i < n; i++) if (fscanf(f, "%d", &c[i]) != 1) exit(1);
    fclose(f); return c;
}
int main(int argc, char** argv) {
    if (const char* s = getenv("K")) {
        char* p = (char*)s;
        while (*p) { char* q; long i = strtol(p, &q, 10); if (q == p || *q != '=') break; p = q + 1; double v = strtod(p, &q); if (q == p) break; if (i >= 0 && i < 64) K[i] = v; p = q; }
    }
#ifdef LWMUT
    if (const char* s = getenv("SFW")) { char* p = (char*)s; int i = 0; while (*p && i < 27) { char* q; double v = strtod(p, &q); if (q == p) break; LW[i++] = v; p = q; while (*p == ',' || *p == ' ') p++; } }
#endif
    const char* dumpAt = getenv("DUMPAT"); long long dumpCi = dumpAt ? atoll(dumpAt) : -1; bool digs = getenv("DIGS") != nullptr;
    int ncut = getenv("SFCUTS") ? atoi(getenv("SFCUTS")) : 0;
    for (int a = 1; a < argc; a++) {
        double w; std::vector<int> cmd = readTest(argv[a], w);
        Alloc* A = new Alloc();
        ORD.assign(1, 0); { for (size_t i = 0; i < cmd.size(); i++) { if (cmd[i] > 0) ORD.push_back(1 << 30); else if (cmd[i] < 0) ORD[-cmd[i]] = (int)i; } }
        FILE* posF = getenv("DUMPPOS") ? fopen(getenv("DUMPPOS"), "w") : nullptr;
        FILE* frontF = getenv("FRONT") ? fopen(getenv("FRONT"), "w") : nullptr; long long frontStep = getenv("FRONTSTEP") ? atoll(getenv("FRONTSTEP")) : 100;
        std::vector<int> cells; std::vector<int> plen(1, 0); std::vector<char> ever; std::vector<int> freedAt; std::vector<int> lastl;
        unsigned long long hsh = 1469598103934665603ULL;   // FNV-хеш всех ответов (сверка «ответы не изменились»)
        ll occ = 0; bool ok = true; double cutSum = 0; int cutCnt = 0;
        long long n = cmd.size();
        for (long long ci = 0; ci < n && ok; ci++) {
            int c = cmd[ci]; CUR = ci; NCMD = (long long)cmd.size();
            if (c == 0) break;
            if (ci == dumpCi) { FILE* f = fopen(getenv("DUMPF"), "w"); for (size_t i = 0; i < cells.size(); i++) fprintf(f, "%d\n", cells[i]); fclose(f); }
            if (frontF && ci % frontStep == 0) fprintf(frontF, "%lld %lld %lld %d\n", ci, (ll)cells.size(), occ, (int)A->tailStart());
            if (ncut > 1 && ci > 0) { long long lo = n - n / 10, st = std::max(1LL, n / 10 / ncut); if (ci >= lo && (ci - lo) % st == 0 && !cells.empty()) { cutSum += (double)occ / cells.size(); cutCnt++; } }
            if (c > 0) {
                ll N, D; A->place(c, N, D);
                if (posF) fprintf(posF, "%lld ", D);
                hsh = (hsh ^ (unsigned long long)(N * 1000003 + D)) * 1099511628211ULL;
                ll T = (ll)cells.size() + 20 * N;
                if (digs && N > 0) printf("dig ci=%lld L=%d N=%lld T=%lld occ=%lld free=%lld\n", ci, c, N, T, occ, (ll)cells.size() - occ);
                if (digs && N > 0) { size_t i = 0; printf("   holes:"); while (i < cells.size()) { if (cells[i]) { i++; continue; } size_t j = i; while (j < cells.size() && !cells[j]) j++; printf(" %zu%s", j - i, j == cells.size() ? "T" : ""); i = j; } printf("\n"); }
                if (N < 0 || T != A->T || D < 0 || D + c > T || T > 100000000) { ok = false; break; }
                cells.resize(T, 0); freedAt.resize(T, -1); ever.resize(T, 0); lastl.resize(T, 0);
                int pid = (int)plen.size(); plen.push_back(c);
                for (int i = 0; i < c; i++) { if (cells[D + i]) { ok = false; break; } cells[D + i] = pid; ever[D + i] = 1; lastl[D + i] = c; }
                occ += c;
            } else {
                int pid = -c; ll d = A->remove(pid);
                hsh = (hsh ^ (unsigned long long)d) * 1099511628211ULL;
                if (d < 0 || d >= (ll)cells.size() || cells[d] != pid) { ok = false; break; }
                for (int i = 0; i < plen[pid]; i++) cells[d + i] = 0, freedAt[d + i] = (int)ci;
                occ -= plen[pid];
            }
        }
        double kpi = ok && !cells.empty() ? (double)occ / cells.size() : 0;
        if (ok && cutCnt > 0) kpi = (cutSum + kpi) / (cutCnt + 1);
        if (getenv("HOLES") && ok) {
            // финальные дыры: размер, возраст (команд с последнего освобождения), длина бывшей трубы
            size_t i = 0, e = cells.size(); while (e > 0 && !cells[e - 1]) e--;
            FILE* f = fopen(getenv("HOLES"), "w");
            while (i < e) { if (cells[i]) { i++; continue; } size_t j = i; int fa = -1; while (j < e && !cells[j]) { if (freedAt[j] > fa) fa = freedAt[j]; j++; }
                fprintf(f, "%zu %lld %d %d\n", j - i, fa < 0 ? -1 : n - fa, lastl[i], (int)ever[i]); i = j; }
            fclose(f);
        }
        if (getenv("DIAG") && ok) {
            long long nev = 0, byz[64] = {0}, byl[64] = {0}, tail = 0; size_t e = cells.size();
            while (e > 0 && !cells[e - 1]) e--;
            tail = cells.size() - e;
            for (size_t i = 0; i < e; i++) if (!cells[i]) { if (!ever[i]) nev++; else byl[lastl[i]]++; int z = A->zl.get(i); byz[z < 63 ? z : 63]++; }
            printf("DIAG tail=%lld never=%lld wasL:", tail, nev);
            for (int l = 1; l <= 40; l++) if (byl[l]) printf(" %d:%lld", l, byl[l]);
            printf(" | zl:"); for (int l = 0; l < 64; l++) if (byz[l]) printf(" %d:%lld", l, byz[l]);
            printf("\n");
            if (getenv("NEVSEG")) { size_t i = 0; printf("NEV:"); while (i < e) { if (cells[i] || ever[i]) { i++; continue; } size_t j = i; while (j < e && !cells[j] && !ever[j]) j++; printf(" %zu+%zu", i, j - i); i = j; } printf(" T=%zu\n", cells.size()); }
        }
        if (frontF) fclose(frontF);
        if (posF) fclose(posF);
        printf("%s w=%g kpi=%.6f T=%zu h=%016llx%s\n", argv[a], w, kpi, cells.size(), hsh, ok ? "" : " WA");
        fflush(stdout);
        delete A;
    }
}
