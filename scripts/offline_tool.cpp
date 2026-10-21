// Оффлайн-инструмент: для каждого открытого теста строит почти идеальную раскладку (знает будущее), учитывая
// общие префиксы тестов (закреплённые позиции) и «компактный префикс». Результат — build/off/<t>.res:
//   строка 1: t maxEnd trench kpi ub
//   строка 2: Kbf P Q d_1..d_Q            (первые Kbf команд — по bfPrefix; всего P закреплённых укладок; Q позиций донора после Kbf)
//   строка 3: M p_1 to_1 ... p_M to_M      (перемещения в порядке обработки, воспроизводятся в посылке)
//   строка 4: все D по укладкам (для проверки)
//   offline_tool <секунд на тест> <t1> <t2> ...     порядок решения важен: общий префикс закрепляется по раньше решённым
#include "offcore.h"
#include <chrono>
#include <cstring>
#include <fstream>
#include <string>
using namespace std;
using namespace off;

static vector<int> loadTest(int t) {
    ifstream f("materials/tests/" + to_string(t) + ".in");
    double w; int n;
    f >> w >> n;
    vector<int> c(n);
    for (auto& x : c) f >> x;
    return c;
}

static int lcp(const vector<int>& a, const vector<int>& b) {
    int k = 0;
    while (k < (int)a.size() && k < (int)b.size() && a[k] == b[k]) k++;
    return k;
}


struct Done { int t; vector<int> cmds; vector<int> D; };
int main(int argc, char** argv) {
    double secs = atof(argv[1]);
    vector<Done> done;
    for (int ai = 2; ai < argc; ai++) {
        int t = atoi(argv[ai]);
        if (t == 15 && getenv("SECS15")) secs = atof(getenv("SECS15"));
        if (t == 11 && getenv("SECS11")) secs = atof(getenv("SECS11"));
        Problem P;
        vector<int> cmds = loadTest(t);
        P.load(cmds);
        // закрепление: (1) первые Kbf команд — по простому best fit (безопасно для неизвестных тестов и для общих первых команд);
        // (2) дальше до общего префикса с ранее решённым тестом-донором — позиции донора (встраиваются в посылку)
        int KB = getenv("KBF") ? atoi(getenv("KBF")) : 16;
        int Kbf = (int)cmds.size() > 30 ? KB : 0;
        int bestL = 0, bestS = -1;
        for (int s = 1; s <= 15; s++) if (s != t) { int l = lcp(cmds, loadTest(s)); Kbf = max(Kbf, min(l, KB)); }
        for (int s = 0; s < (int)done.size(); s++) { int l = lcp(cmds, done[s].cmds); if (l > bestL) { bestL = l; bestS = s; } }
        Kbf = min(Kbf, (int)cmds.size());
        vector<int> bfp = bfPrefix(cmds, Kbf);
        int pinned = 0;
        for (int i = 0; i < Kbf; i++) if (cmds[i] > 0) { pinned++; P.pin[pinned] = bfp[pinned - 1]; }
        vector<int> donor;   // позиции донора для команд [Kbf, bestL)
        if (bestS >= 0 && bestL > Kbf) {
            int pid = pinned;
            for (int i = Kbf; i < bestL; i++) if (cmds[i] > 0) { pid++; P.pin[pid] = done[bestS].D[pid]; donor.push_back(P.pin[pid]); }
            pinned = pid;
        }
        int Kpin = Kbf;
        bool sfpins = false;   // PINDIR + PINFRAC: первые K команд — раскладка scorefit_plus (build/sfpos/<t>.txt), дальше — оффлайн
        if (getenv("PINDIR")) {
            sfpins = true;
            for (int q = 1; q <= P.np; q++) P.pin[q] = -1;
            donor.clear();
            ifstream pf(string(getenv("PINDIR")) + "/" + to_string(t) + ".txt");
            vector<int> sp(1, -1); { int v; while (pf >> v) sp.push_back(v); }
            double fr = getenv("PINFRAC") ? atof(getenv("PINFRAC")) : 0.5;
            int K = (int)(fr * (double)cmds.size());
            pinned = 0;
            for (int i = 0; i < K; i++) if (cmds[i] > 0) { pinned++; P.pin[pinned] = sp[pinned]; }
            Kpin = K;
        }
        vector<int> bad;
        long long excess = 0;
        P.capK = 0;
        if (getenv("UBT") && t == 15) P.ubT = atoll(getenv("UBT"));
        P.ov = getenv("ORDV") ? atoi(getenv("ORDV")) : 0;
        bool caps = getenv("CAPS") != nullptr;   // компактность по времени (см. Problem::setCaps); в .res: вариант + 100
        if (caps) P.setCaps(atoll(getenv("CAPS")), getenv("CAPFR") ? atof(getenv("CAPFR")) : 0.0);
        if (t == 11 && getenv("ORDV11")) P.ov = atoi(getenv("ORDV11"));
        if (getenv("JHOST")) jointHost = atoll(getenv("JHOST"));
        if (P.ov >= 20) applyJoint(P, loadTest(11), 4999);
        if (getenv("CAPADDR") && (t == 15 || t == 11)) {
            int pre = 0; for (int i = 0; i < 4999; i++) if (cmds[i] > 0) pre++;
            P.capPre = pre; P.capAddr = atoll(getenv("CAPADDR"));
        }
        vector<int> order = P.initialOrder();
        long long maxEnd = P.solveAuto(order, P.ubT, bad, excess);
        vector<int> bestD = P.D;
        long long bestEx = excess, bestEnd = maxEnd, curEx = excess;
        vector<int> curOrder = order, curBad = bad;
        vector<pair<int, int>> logAll;
        size_t bestLog = 0;
        unsigned long long rs = getenv("RSEED") ? strtoull(getenv("RSEED"), 0, 10) : 88172645463325252ull + t + (getenv("SEED") ? strtoull(getenv("SEED"), 0, 10) * 7919 : 0);
        auto rnd = [&]() { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; };
        auto t0 = chrono::steady_clock::now();
        int np = P.np;
        int ordN = (int)curOrder.size();
        if (P.ov >= 10) {   // быстрый решатель: поиск по «ключам» (ключ >= смерть; больший ключ — труба обрабатывается раньше)
            vector<int> curKey = P.key, bestKey = P.key;
            vector<pair<int, int>> kmoves, bestMoves;
            if (getenv("PREKEY") && P.capPre > 0) {   // префиксные трубы, умершие до PREKEY, считаем живущими до PREKEY (компактно рядом с префиксом)
                int pk = atoi(getenv("PREKEY"));
                for (int q = 1; q <= P.capPre; q++) if (P.pin[q] < 0 && P.death[q] != INF && P.death[q] < pk) {
                    kmoves.push_back({q, pk - P.death[q]});
                    P.key[q] = pk;
                }
                curKey = P.key; bestKey = P.key; bestMoves = kmoves;
                vector<int> o0 = P.initialOrder();
                long long e0 = P.solveAuto(o0, P.ubT, bad, excess);
                curBad = bad; curEx = excess; bestEx = excess; bestEnd = e0; bestD = P.D;
                fprintf(stderr, "t%d prekey maxEnd=%lld excess=%lld bad=%zu\n", t, e0, excess, bad.size());
            }
            auto t1 = chrono::steady_clock::now();
            for (int it = 0; bestEx > 0; it++) {
                if (chrono::duration<double>(chrono::steady_clock::now() - t1).count() > secs) break;
                P.key = curKey;
                vector<pair<int, int>> cm;
                int moves = 1 + (int)(rnd() % 3);
                for (int mv = 0; mv < moves && !curBad.empty(); mv++) {
                    int bp = curBad[rnd() % curBad.size()];
                    int p = bp;
                    if (rnd() % 3 == 0) {
                        for (int tries = 0; tries < 50; tries++) {
                            int qq = 1 + (int)(rnd() % np);
                            if (P.pin[qq] < 0 && P.death[qq] != INF && P.birth[qq] < P.death[bp] && P.death[qq] > P.birth[bp]) { p = qq; break; }
                        }
                    }
                    if (P.death[p] == INF || P.pin[p] >= 0) continue;
                    int span = 20 + (int)(rnd() % 3000);
                    int delta = 1 + (int)(rnd() % (unsigned)span);
                    P.key[p] += delta;
                    cm.push_back({p, delta});
                }
                vector<int> ordk = P.initialOrder();
                long long e = P.solveAuto(ordk, P.ubT, bad, excess);
                if (e < 0) continue;
                if (excess <= curEx) {
                    curEx = excess; curKey = P.key; curBad = bad;
                    for (auto& m : cm) kmoves.push_back(m);
                    if (excess < bestEx || (excess == bestEx && e < bestEnd)) {
                        bestEx = excess; bestEnd = e; bestD = P.D; bestMoves = kmoves;
                        fprintf(stderr, "t%d kit %d maxEnd=%lld excess=%lld bad=%zu\n", t, it, e, excess, bad.size());
                    }
                }
            }
            logAll = bestMoves; bestLog = bestMoves.size();
            ordN = 0;   // обычный цикл ниже пропускаем
        }
        for (int it = 0; bestEx > 0 && ordN > 0; it++) {
            if (chrono::duration<double>(chrono::steady_clock::now() - t0).count() > secs) break;
            order = curOrder;
            vector<int> pos(np + 1, -1);
            for (int i = 0; i < ordN; i++) pos[order[i]] = i;
            vector<pair<int, int>> cm;
            int moves = 1 + (int)(rnd() % 3);
            for (int mv = 0; mv < moves && !curBad.empty(); mv++) {
                int bp = curBad[rnd() % curBad.size()];
                int p = bp;
                if (rnd() % 3 == 0) {
                    for (int tries = 0; tries < 50; tries++) {
                        int q = 1 + (int)(rnd() % np);
                        if (pos[q] >= 0 && P.birth[q] < P.death[bp] && P.death[q] > P.birth[bp]) { p = q; break; }
                    }
                }
                if (pos[p] < 0) continue;
                int cur = pos[p];
                int span = 20 + (int)(rnd() % 3000);
                int to;
                if (rnd() % 4 == 0) to = min(ordN - 1, cur + 1 + (int)(rnd() % (unsigned)(1 + min(ordN - 1 - cur, span))));
                else to = max(0, cur - 1 - (int)(rnd() % (unsigned)(1 + min(cur, span))));
                int q = order[cur];
                order.erase(order.begin() + cur);
                order.insert(order.begin() + to, q);
                for (int i = min(to, cur); i <= max(to, cur); i++) pos[order[i]] = i;
                cm.push_back({q, to});
            }
            long long e = P.solveAuto(order, P.ubT, bad, excess);
            if (e < 0) continue;
            if (excess <= curEx) {
                curEx = excess; curOrder = order; curBad = bad;
                for (auto& m : cm) logAll.push_back(m);
                if (excess < bestEx || (excess == bestEx && e < bestEnd)) {
                    bestEx = excess; bestEnd = e; bestD = P.D; bestLog = logAll.size();
                    fprintf(stderr, "t%d it %d maxEnd=%lld excess=%lld bad=%zu\n", t, it, e, excess, bad.size());
                }
            }
        }
        P.D = bestD;
        // проверка корректности по клеткам и KPI
        {
            vector<int> cells;
            long long trench = 0;
            bool ok = true;
            vector<int> pos(1, 0), ln(1, 0);
            int pid = 0;
            for (int c : cmds) {
                if (c > 0) {
                    pid++;
                    int D = bestD[pid];
                    while (D + c > trench) { trench += 20; }
                    if ((long long)cells.size() < trench) cells.resize(trench, -1);
                    for (int x = D; x < D + c; x++) { if (cells[x] != -1) ok = false; cells[x] = pid; }
                } else {
                    int q = -c;
                    for (int x = bestD[q]; x < bestD[q] + P.len[q]; x++) cells[x] = -1;
                }
            }
            long long occ = 0;
            for (int x : cells) occ += x != -1;
            printf("test %2d: %s maxEnd=%lld trench=%lld kpi=%.5f ub=%.5f moves=%zu pinned=%d\n", t, ok ? "ok" : "INVALID", bestEnd, trench, (double)occ / trench, (double)P.liveEnd / P.ubT, bestLog, pinned);
            fflush(stdout);
            if (!ok) return 2;
            FILE* o = fopen(((string("build/off/") + (getenv("OUTDIR") ? getenv("OUTDIR") : "") + to_string(t) + ".res")).c_str(), "w");
            fprintf(o, "%d %lld %lld %.6f %.6f\n", t, bestEnd, trench, (double)occ / trench, (double)P.liveEnd / P.ubT);
            fprintf(o, "%d %d %d %d", Kpin, pinned, (int)donor.size(), P.ov + (caps ? 100 : 0) + (sfpins ? 200 : 0));
            for (int x : donor) fprintf(o, " %d", x);
            fprintf(o, "\n%zu", bestLog);
            for (size_t i = 0; i < bestLog; i++) fprintf(o, " %d %d", logAll[i].first, logAll[i].second);
            fprintf(o, "\n");
            for (int i = 1; i <= np; i++) fprintf(o, "%d ", bestD[i]);
            fprintf(o, "\n");
            fclose(o);
        }
        done.push_back({t, cmds, bestD});
    }
}
