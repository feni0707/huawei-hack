// scorefit: размещение трубы по функции стоимости (см. docs/variants/scorefit.md).
// Сборка для посылки: c++ -std=c++17 -O2 -o build/bin/scorefit solutions/scorefit.cpp
// Сборка оффлайн-симулятора (без интерактора, для экспериментов):
//   c++ -std=c++17 -O2 -DSIM -o build/sim_scorefit solutions/scorefit.cpp
//   build/sim_scorefit materials/tests/*.in      (параметры: env SF="p0 p1 ...")
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <unistd.h>

typedef long long ll;

// ---------------------------------------------------------------- параметры
// P0 затухание статистики длин, P1 порог доли для minRecent, P2 цена мёртвого метра остатка,
// P3 цена живого метра остатка, P4 сколько промежутков каждой длины смотреть, P5 бонус точного попадания,
// P6 цена метра остатка (best fit), P7 порог стоимости, выше которого кладём в хвост (0 — никогда),
// P8 выбор стороны (0 лево, 1 остаток к молодому соседу, 2 к старому, 3 право)
// P10 потоки (0 нет, 1 оракул, 2 по признакам), P11 горизонт оракула, P20 размер участка для потока B, м,
// P21 мин. число пар удалений для включения, P22 порог каппы, P23 затухание рождений, P24 затухание пар
static double P[32] = {0.999, 0.002, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       200, 20, 0.6, 0.999, 0.998};

// ---------------------------------------------------------------- аллокатор
// Свободные промежутки: граничные метки (склейка за O(1)) + корзины по точной длине (< NB) + set для больших.
static const int NB = 128;
static const int MAXL = 40;

struct Alloc {
    ll T = 0;                          // длина траншеи
    std::vector<int> gapLenAt;         // длина промежутка, начинающегося в s (0 — нет)
    std::vector<int> gapStartAtEnd;    // начало промежутка, кончающегося в e (-1 — нет)
    std::vector<int> pidStartAt;       // PID трубы, начинающейся в x (0 — нет)
    std::vector<int> pidEndAt;         // PID трубы, кончающейся в x (0 — нет)
    std::set<int> bucket[NB];          // начала промежутков длины s (s < NB), по адресу
    std::set<std::pair<int,int>> big;  // (len, start) для len >= NB
    std::vector<int> pos, len;         // по PID
    int nPlaced = 0;

    // онлайн-статистика длин (с затуханием)
    double cnt[MAXL + 1] = {0}, cntTot = 0;

    void ensure(ll n) {
        if ((ll)gapLenAt.size() < n + 1) {
            ll m = std::max<ll>(n + 1, gapLenAt.size() * 2 + 64);
            gapLenAt.resize(m, 0); gapStartAtEnd.resize(m, -1);
            pidStartAt.resize(m, 0); pidEndAt.resize(m, 0);
        }
    }
    void addGap(int s, int l) {
        gapLenAt[s] = l; gapStartAtEnd[s + l] = s;
        if (l < NB) bucket[l].insert(s); else big.insert({l, s});
    }
    void delGap(int s) {
        int l = gapLenAt[s];
        if (l < NB) bucket[l].erase(s); else big.erase({l, s});
        gapLenAt[s] = 0; gapStartAtEnd[s + l] = -1;
    }
    int freeRange(int s, int l) {   // возвращает конец склеенного промежутка
        int e = s + l;
        if (e < T && gapLenAt[e] > 0) { int r = gapLenAt[e]; delGap(e); e += r; }
        if (gapStartAtEnd[s] >= 0) { int ps = gapStartAtEnd[s]; delGap(ps); s = ps; }
        addGap(s, e - s);
        return e;
    }
    void takeFrom(int gs, int d, int l) {
        int gl = gapLenAt[gs];
        delGap(gs);
        if (d > gs) addGap(gs, d - gs);
        if (d + l < gs + gl) addGap(d + l, gs + gl - (d + l));
    }
    void dig(ll n) {
        if (n <= 0) return;
        ll oldT = T;
        T += 20 * n;
        ensure(T);
        freeRange((int)oldT, (int)(T - oldT));
    }
    int tailStart() const { return (T > 0 && gapStartAtEnd[T] >= 0) ? gapStartAtEnd[T] : (int)T; }

    // ---- статистика
    void observe(int L) {
        double dec = P[0];
        for (int l = 1; l <= MAXL; l++) cnt[l] *= dec;
        cntTot = cntTot * dec + 1; cnt[L] += 1;
    }
    // минимальная длина, которую «ещё ждём» (доля в недавних укладках > eps)
    int minRecent() const {
        double eps = P[1];
        for (int l = 1; l <= MAXL; l++) if (cnt[l] > eps * cntTot) return l;
        return 1;
    }

    // стоимость остатка r после укладки (r > 0)
    double fragCost(int r, int m) const {
        if (r < m) return r * P[2];        // мёртвый остаток
        return r * P[3];                                           // живой остаток (P3=0 — без разницы)
    }

    // ---- потоки (сегрегация по признаку, чьи носители умирают вместе)
    // stream[pid]: -1 — обычная, 0 — поток A (кладётся к левому краю граничного промежутка),
    // 1 — поток B (к правой границе роста F, растёт влево)
    std::vector<signed char> stream;
    int F = 0;                        // правая граница роста; F == T — обычный режим
    std::vector<unsigned char> feat;  // биты признаков при рождении: bit0 — чётность индекса в серии равных длин, bit1 — чётность PID
    int runIdx = 0, lastL = 0;
    static const int NF = 2;
    double coSame[NF] = {0}, coTot[NF] = {0};   // соседние удаления в пачке: совпал ли признак
    double birthV[NF][2] = {{0}}, deathV[NF][2] = {{0}};
    int lastDelPid = 0; bool lastWasDel = false;
    const int* oracleDeath = nullptr;   // только SIM-эксперименты: время смерти по PID (команда), -1 — не умрёт
    int curCmd = 0;

    double kappa(int f) const {
        if (coTot[f] < P[21]) return 0;
        double b0 = birthV[f][0], b1 = birthV[f][1], bt = b0 + b1;
        if (bt <= 0) return 0;
        double e = (b0 * b0 + b1 * b1) / (bt * bt);   // ожидаемая доля совпадений при независимости
        double o = coSame[f] / coTot[f];
        return e >= 1 ? 0 : (o - e) / (1 - e);
    }
    double pdie(int f, int v) const { return (deathV[f][v] + 1) / (birthV[f][v] + 2); }
    int chooseStream(int pid, int L) {
        unsigned char fv = feat[pid];
        if (P[10] == 1 && oracleDeath) {
            int dt = oracleDeath[pid];
            return dt >= 0 && dt - curCmd <= P[11] ? 1 : 0;
        }
        if (P[10] != 2) return -1;
        // признак: каппа «умирают вместе» выше порога; при равенстве — тот, где сильнее разница p_die
        int bf = -1; double bs = -1;
        for (int f = 0; f < NF; f++) {
            double k = kappa(f);
            if (k < P[22]) continue;
            double s = k + std::fabs(pdie(f, 0) - pdie(f, 1));
            if (s > bs) { bs = s; bf = f; }
        }
        if (bf < 0) return -1;
        int v = (fv >> bf) & 1;
        // тот, кто умирает чаще, — в поток B (правый край, у границы роста)
        double pv = pdie(bf, v), po = pdie(bf, v ^ 1);
        if (std::fabs(pv - po) < 0.1) return v;
        return pv > po ? 1 : 0;
    }

    // положить в промежуток gs длины s: выбрать сторону, вернуть D
    int chooseSide(int gs, int s, int L, int st) {
        int r = s - L;
        if (r == 0) return gs;
        int lp = pidEndAt[gs], rp = pidStartAt[gs + s];
        if (st >= 0) {
            // прижимаем к соседу из того же потока
            int lc = lp ? stream[lp] : -2, rc = rp ? stream[rp] : -2;
            if (lc == st && rc != st) return gs;
            if (rc == st && lc != st) return gs + r;
            return st == 1 ? gs + r : gs;
        }
        if (P[8] == 1) { if (lp > rp) return gs + r; }
        else if (P[8] == 3) return gs + r;
        return gs;
    }

    std::pair<ll,ll> place(int L) {
        int m = minRecent();
        int pid = (int)pos.size();
        if (L == lastL) runIdx++; else runIdx = 0;
        lastL = L;
        feat.push_back((unsigned char)((runIdx & 1) | ((pid & 1) << 1)));
        int st = chooseStream(pid, L);
        stream.push_back((signed char)st);
#ifdef SIM
        if (getenv("SFDBG") && pid % atoi(getenv("SFDBG")) == 0)
            fprintf(stderr, "pid=%d L=%d feat=%d st=%d k0=%.3f k1=%.3f tot=%.1f,%.1f F=%d T=%lld\n", pid, L, feat[pid], st, kappa(0), kappa(1), coTot[0], coTot[1], F, T);
#endif
        for (int f = 0; f < NF; f++) {
            birthV[f][0] *= P[23]; birthV[f][1] *= P[23];
            birthV[f][(feat[pid] >> f) & 1] += 1;
        }
        if (F > T) F = (int)T;
        int ts = tailStart();
        double bestC = 1e18; int bgs = -1, bs = 0;
        int M = P[4] > 0 ? (int)P[4] : 1;  // сколько промежутков каждой длины смотреть
        auto eval = [&](int gs, int s) {
            double cst;
            if (P[14] == 1) {
                cst = gs;                                   // first fit по адресу
            } else if (P[14] == 2) {
                cst = -s;                                   // worst fit
            } else if (P[14] == 3) {
                cst = s == L ? -1e9 : -s;                   // точное, иначе worst
            } else if (P[14] == 4) {
                cst = -(double)gs;                          // last fit
            } else {
                cst = s == L ? -P[5] : fragCost(s - L, m) + P[6] * (s - L);
            }
            if (cst < bestC) { bestC = cst; bgs = gs; bs = s; }
        };
        for (int s = L; s < NB; s++) {
            int k = 0;
            for (auto it = bucket[s].begin(); it != bucket[s].end() && k < M; ++it) {
                if (*it == ts || *it + s == F) continue;  // хвост и граничный промежуток — отдельно
                k++; eval(*it, s);
            }
        }
        {
            int k = 0;
            for (auto it = big.begin(); it != big.end() && k < M; ++it) {
                if (it->second == ts || it->second + it->first == F) continue;
                k++; eval(it->second, it->first);
            }
        }
        ll N = 0; int D;
        bool useGap = bgs >= 0 && bestC < P[7] + 1e17 * (P[7] == 0);
        if (st == 1 && P[12] == 2) useGap = false;
        if (useGap) {
            D = chooseSide(bgs, bs, L, st);
            takeFrom(bgs, D, L);
        } else {
            int fgs = (F > 0 && gapStartAtEnd[F] >= 0) ? gapStartAtEnd[F] : F;
            int fl = F - fgs;
            int chunk = (int)P[20];
            if (fl >= L) {
                D = st == 1 ? F - L : fgs;
                takeFrom(fgs, D, L);
                if (st == 1) F = D;
            } else if (F == T) {
                int need = L - fl;
                if (st == 1) need = std::max(need, chunk - fl);
                N = (need + 19) / 20;
                dig(N);
                D = st == 1 ? (int)T - L : fgs;
                takeFrom(fgs, D, L);
                F = st == 1 ? D : (int)T;
            } else {
                // граничный промежуток исчерпан — новый участок в конце траншеи
                ts = tailStart();
                int tl = (int)(T - ts);
                int need = L - tl;
                if (st >= 0) need = std::max(need, chunk - tl);
                N = need > 0 ? (need + 19) / 20 : 0;
                dig(N);
                ts = tailStart();
                D = st == 1 ? (int)T - L : ts;
                takeFrom(ts, D, L);
                F = st == 1 ? D : (int)T;
            }
        }
        pidStartAt[D] = pid; pidEndAt[D + L] = pid;
        pos.push_back(D); len.push_back(L);
        observe(L); nPlaced++;
        lastWasDel = false;
        return {N, D};
    }
    ll remove(int pid) {
        int d = pos[pid], l = len[pid];
        pidStartAt[d] = 0; pidEndAt[d + l] = 0;
        int e = freeRange(d, l);
        if (d == F) F = e;
        unsigned char fv = feat[pid];
        for (int f = 0; f < NF; f++) {
            int v = (fv >> f) & 1;
            deathV[f][v] += 1;
            if (lastWasDel) {
                coSame[f] *= P[24]; coTot[f] *= P[24];
                coTot[f] += 1; coSame[f] += (v == ((feat[lastDelPid] >> f) & 1));
            }
        }
        lastWasDel = true; lastDelPid = pid;
        return d;
    }
    size_t gapCount() const { size_t c = big.size(); for (int s = 0; s < NB; s++) c += bucket[s].size(); return c; }
    void holeStats(ll& holes, ll& small) const {
        holes = small = 0; int ts = tailStart();
        for (int s = 1; s < NB; s++) for (int g : bucket[s]) if (g != ts) { holes += s; if (s < 8) small += s; }
        for (auto& g : big) if (g.second != ts) holes += g.first;
    }
    Alloc() { pos.push_back(-1); len.push_back(0); feat.push_back(0); stream.push_back(-1); ensure(1024); }
};

#ifdef SIM
// ---------------------------------------------------------------- оффлайн-симулятор
// Читает тест целиком (это только для экспериментов: аллокатор видит команды строго по одной),
// проверяет корректность как интерактор и печатает KPI.
static std::vector<int> readTest(const char* fn, double& w) {
    FILE* f = fopen(fn, "r");
    if (!f) { perror(fn); exit(1); }
    int n; if (fscanf(f, "%lf %d", &w, &n) != 2) exit(1);
    std::vector<int> c(n);
    for (int i = 0; i < n; i++) if (fscanf(f, "%d", &c[i]) != 1) exit(1);
    fclose(f);
    return c;
}
int main(int argc, char** argv) {
    if (const char* s = getenv("SF")) {
        // формат: "i=v i=v ..."
        char* p = (char*)s;
        while (*p) {
            char* q; long i = strtol(p, &q, 10); if (q == p || *q != '=') break;
            p = q + 1; double v = strtod(p, &q); if (q == p) break; if (i >= 0 && i < 32) P[i] = v; p = q;
        }
    }
    bool quiet = getenv("SFQ") != nullptr;
    double total = 0, wsum = 0;
    for (int a = 1; a < argc; a++) {
        double w; std::vector<int> cmd = readTest(argv[a], w);
        Alloc A;
        std::vector<int> death(1, -1);
        { int p = 0; for (size_t i = 0; i < cmd.size(); i++) { if (cmd[i] > 0) { death.push_back(-1); p++; } else death[-cmd[i]] = (int)i; } }
        A.oracleDeath = death.data();
        std::vector<int> cells;  // pid или 0
        std::vector<int> lastLen;
        ll occ = 0; bool ok = true;
        long long step = getenv("SFTRACE") ? atoll(getenv("SFTRACE")) : 0; long long ci = 0;
        for (int c : cmd) {
            if (step && ci % step == 0) {
                // гистограмма дыр
                long long holes, small; A.holeStats(holes, small);
                printf("  t=%lld T=%lld occ=%lld waste=%lld holes=%lld small<8=%lld\n", ci, A.T, occ, A.T - occ, holes, small);
            }
            A.curCmd = (int)ci;
            ci++;
            if (c > 0) {
                auto r = A.place(c);
                ll T = (ll)cells.size() + 20 * r.first;
                if (r.first < 0 || T != A.T || r.second < 0 || r.second + c > T) { ok = false; break; }
                cells.resize(T, 0); lastLen.resize(T, 0);
                for (int i = 0; i < c; i++) if (cells[r.second + i]) { ok = false; break; }
                if (!ok) break;
                int pid = (int)A.pos.size() - 1;
                for (int i = 0; i < c; i++) cells[r.second + i] = pid, lastLen[r.second + i] = c;
                occ += c;
            } else {
                int pid = -c;
                ll d = A.remove(pid);
                if (d < 0 || d >= (ll)cells.size() || cells[d] != pid) { ok = false; break; }
                int l = A.len[pid];
                for (int i = 0; i < l; i++) cells[d + i] = 0;
                occ -= l;
            }
        }
        if (getenv("SFORIG") && ok) {
            // происхождение свободных клеток в конце: 0 — не занималась, иначе длина последней трубы на этом месте
            long long byLen[41] = {0}, never = 0;
            for (size_t i = 0; i < cells.size(); i++) if (!cells[i]) { if (lastLen[i]) byLen[lastLen[i]]++; else never++; }
            printf("  free: never=%lld", never);
            long long sh = 0, lg = 0; for (int l = 1; l <= 40; l++) (l <= 7 ? sh : lg) += byLen[l];
            printf(" was_short=%lld was_long=%lld\n", sh, lg);
        }
        double kpi = ok && !cells.empty() ? (double)occ / cells.size() : 0;
        if (!ok) fprintf(stderr, "WA on %s\n", argv[a]);
        total += kpi * w; wsum += w;
        if (!quiet) printf("%-40s w=%-4g kpi=%.5f T=%zu occ=%lld\n", argv[a], w, kpi, cells.size(), occ);
    }
    printf("TOTAL %.5f / %.4g  avg %.5f\n", total, wsum, total / wsum);
}
#else
// ---------------------------------------------------------------- интерактивный режим
static char ibuf[1 << 16];
static int ipos = 0, ilen = 0;
static inline int gc() {
    if (ipos == ilen) { ilen = (int)read(0, ibuf, sizeof ibuf); ipos = 0; if (ilen <= 0) return -1; }
    return ibuf[ipos++];
}
static bool readInt(int& x) {
    int c = gc();
    while (c != '-' && (c < '0' || c > '9')) { if (c == -1) return false; c = gc(); }
    bool neg = false; if (c == '-') { neg = true; c = gc(); }
    x = 0;
    while (c >= '0' && c <= '9') { x = x * 10 + (c - '0'); c = gc(); }
    if (neg) x = -x;
    return true;
}
static void writeOut(const char* s, int n) {
    while (n > 0) { ssize_t k = write(1, s, n); if (k <= 0) return; s += k; n -= (int)k; }
}
static int putNum(char* b, ll v) {
    char t[24]; int n = 0;
    if (v == 0) t[n++] = '0';
    while (v > 0) { t[n++] = char('0' + v % 10); v /= 10; }
    for (int i = 0; i < n; i++) b[i] = t[n - 1 - i];
    return n;
}
int main() {
    Alloc A;
    int c; char ob[64];
    while (readInt(c)) {
        if (c == 0) return 0;
        int n;
        if (c > 0) {
            auto r = A.place(c);
            n = putNum(ob, r.first); ob[n++] = ' '; n += putNum(ob + n, r.second);
        } else {
            n = putNum(ob, A.remove(-c));
        }
        ob[n++] = '\n';
        writeOut(ob, n);
    }
    return 0;
}
#endif
