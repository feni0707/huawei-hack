// scorefit_lab: ЭКСПЕРИМЕНТАЛЬНАЯ версия scorefit со всеми ручками P[] и оффлайн-симулятором (-DSIM).
// Посылать solutions/scorefit.cpp (чистая версия, ответы побитово совпадают с этой при параметрах по умолчанию).
// Сборка оффлайн-симулятора (без интерактора, для экспериментов):
//   c++ -std=c++17 -O2 -DSIM -o build/sim_sf solutions/scorefit_lab.cpp
//   build/sim_sf materials/tests/*.in      (параметры: env SF="i=v i=v ...")
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
// P15 мёртвый остаток >= P15 — в хвост (0 выкл), P16 штраф мёртвого остатка с учётом риска смерти соседей
//      0  1      2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
// P13 «крошечная траншея» (T < P13): точное попадание, иначе самая большая дыра; P32 — штраф точного
//      попадания у соседа с высокой долей смертей его длины (только в этом режиме)
// P33=1 — правило P15 только вне потоков, P33=2 — промежутки с мёртвым остатком ≥ P15 исключаются из кандидатов (вне потоков); P36 — мин. доля редкого значения признака для каппы;
// P39 — штраф замуровывания, только если s + len(соседа) встречается среди недавних укладок с долей ≥ P39
static double P[40] = {0.999, 0.002, 0, 0, 1, 0, 1, 0, 1, 0.5, 2, 0, 3, 60, 0, 4, 0, 3, 4, 0,
                       400, 20, 0.6, 0.999, 0.998, 0, 0, 0.99, 0.5, 0, 4, 0.98, 2000, 2, 0, 0, 0.3, 0, 0.99, 0.05};

// ---------------------------------------------------------------- аллокатор
// Свободные промежутки: граничные метки (склейка за O(1)) + корзины по точной длине (< NB) + set для больших.
static const int NB = 128;
static const int MAXL = 40;

struct Alloc {
    ll T = 0;                          // длина траншеи
    // tag[x] > 0 — в x начинается промежуток длины tag[x]; tag[x] < 0 — в x кончается промежуток,
    // начавшийся в -tag[x]-1; 0 — ни то ни другое (соседние промежутки всегда склеены, поэтому одно x
    // не бывает сразу началом и концом). owner[x] — PID трубы в клетке x (в свободных клетках — мусор).
    // Постраничное хранение: память растёт вместе с траншеей без копирования (худший случай T ≈ 8·10^6).
    struct Paged {
        static const int SH = 16, MASK = (1 << SH) - 1;
        std::vector<int*> pg; int fill;
        explicit Paged(int f) : fill(f) {}
        ~Paged() { for (int* p : pg) delete[] p; }
        void grow(ll n) { while ((ll)pg.size() << SH < n) { int* p = new int[1 << SH]; for (int i = 0; i <= MASK; i++) p[i] = fill; pg.push_back(p); } }
        int& operator[](ll i) { return pg[i >> SH][i & MASK]; }
        int operator[](ll i) const { return pg[i >> SH][i & MASK]; }
    };
    Paged tag{0}, owner{0};
    std::set<int> bucket[NB];          // начала промежутков длины s (s < NB), по адресу
    std::set<std::pair<int,int>> big;  // (len, start) для len >= NB
    std::vector<int> pos, len;         // по PID
    int nPlaced = 0;

    // онлайн-статистика длин (с затуханием)
    double cnt[MAXL + 1] = {0}, cntTot = 0;

    void ensure(ll n) {
        tag.grow(n + 1); owner.grow(n + 1);
    }
    void addGap(int s, int l) {
        tag[s] = l; tag[s + l] = -s - 1;
        if (l < NB) bucket[l].insert(s); else big.insert({l, s});
    }
    void delGap(int s) {
        int l = tag[s];
        if (l < NB) bucket[l].erase(s); else big.erase({l, s});
        tag[s] = 0; tag[s + l] = 0;
    }
    int freeRange(int s, int l) {   // возвращает конец склеенного промежутка
        int e = s + l;
        if (e < T && tag[e] > 0) { int r = tag[e]; delGap(e); e += r; }
        if (tag[s] < 0) { int ps = -tag[s] - 1; delGap(ps); s = ps; }
        addGap(s, e - s);
        return e;
    }
    void takeFrom(int gs, int d, int l) {
        int gl = tag[gs];
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
    int tailStart() const { return (T > 0 && tag[T] < 0) ? -tag[T] - 1 : (int)T; }
    int pidLeftOf(int x) const { return x > 0 ? owner[x - 1] : 0; }      // труба, кончающаяся в x (x — начало дыры)
    int pidRightOf(int x) const { return x < T ? owner[x] : 0; }         // труба, начинающаяся в x (x — конец дыры)

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

    // ---- модель ожидаемых потерь дыры: W[r] — сколько метров дыры длины r останется пустыми,
    // если её заполнять случайными трубами из текущего распределения длин
    double W[NB] = {0};
    int wAge = 1 << 30;
    void rebuildW() {
        double pr[MAXL + 1], tot = 0;
        for (int l = 1; l <= MAXL; l++) { pr[l] = cnt[l] > P[1] * cntTot ? cnt[l] : 0; tot += pr[l]; }
        W[0] = 0;
        for (int r = 1; r < NB; r++) {
            double s = 0, w = 0;
            for (int l = 1; l <= MAXL && l <= r; l++) if (pr[l] > 0) { s += pr[l]; w += pr[l] * W[r - l]; }
            W[r] = s > 0 ? w / s : r;
        }
        wAge = 0;
    }
    double Wf(int r) const { return r < NB ? W[r] : 0; }

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
    double bVol = 0;                  // затухающий объём укладок потока B (для размера участка)
    double fbVol = 0;                 // то же, но только укладки B к границе роста
    std::vector<unsigned char> feat;  // биты признаков при рождении: bit0 — чётность индекса в серии равных длин, bit1 — чётность PID
    int runIdx = 0, lastL = 0, curRun = 0;
    std::vector<int> runOf;           // номер серии равных длин по PID
    static const int NF = 2;
    double coSame[NF] = {0}, coTot[NF] = {0};   // соседние удаления в пачке: совпал ли признак
    double birthV[NF][2] = {{0}}, deathV[NF][2] = {{0}};
    int lastDelPid = 0; bool lastWasDel = false;
    int curFeat = -1;
    double tieAge = 0;
    double bornLab[2] = {0, 0}, diedLab[2] = {0, 0};   // проверка прогноза: смерти среди помеченных A/B
    // ---- риск смерти соседа по возрасту
    // aliveBit — дерево Фенвика по PID (живые), deathAge[b] — число смертей в возрасте [2^b-1, 2^(b+1)-1)
    std::vector<int> aliveBit;
    double deathAge[24] = {0};
    void bitAdd(int i, int v) { for (; i < (int)aliveBit.size(); i += i & -i) aliveBit[i] += v; }
    int bitSum(int i) const { int s = 0; for (; i > 0; i -= i & -i) s += aliveBit[i]; return s; }
    static int ageBucket(int a) { int b = 0; while ((2 << b) - 1 <= a && b < 23) b++; return b; }
    // P(труба PID p когда-нибудь умрёт | дожила до текущего возраста)
    double hazard(int p) const {
        if (p <= 0) return 0;   // край траншеи
        int cur = (int)pos.size() - 1;
        int a = cur - p, b = ageBucket(a);
        double dGE = 0; for (int k = b; k < 24; k++) dGE += deathAge[k];
        int lo = cur - ((1 << b) - 1);   // живые с возрастом >= 2^b-1
        double aGE = lo > 0 ? bitSum(lo) : 0;
        return dGE + aGE > 0 ? dGE / (dGE + aGE) : 0.5;
    }
    double birthL[MAXL + 1] = {0}, deathL[MAXL + 1] = {0};
    double pdieLen(int p) const { if (p <= 0) return 0; int l = len[p]; return (deathL[l] + 0.5) / (birthL[l] + 1); }
    // интенсивность смертей по длине: доля недавних удалений длины l / доля живых длины l
    double recDel[MAXL + 1] = {0}, recDelTot = 0; int aliveCntL[MAXL + 1] = {0}, aliveTot = 0;
    double hzRate(int p) const {
        if (p <= 0 || recDelTot < 1) return 0;
        int l = len[p];
        double rel = (recDel[l] / recDelTot) / ((double)std::max(1, aliveCntL[l]) / std::max(1, aliveTot));
        return std::min(1.0, rel / P[30]);
    }
    double hz2(int p) const { return P[18] == 4 ? hzRate(p) : P[18] == 1 ? pdieLen(p) : P[18] == 2 ? std::max(pdieLen(p), hazard(p)) : P[18] == 3 ? pdieLen(p) * hazard(p) : hazard(p); }
#ifdef SIM
    unsigned long long rs = 88172645463325252ULL; int epsRand = 0;
    int forcePid = -1, forceAlt = 0, listPid = -1, optBest = 0; std::vector<int> optList;
    unsigned long long rng() { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
#endif
    const int* oracleDeath = nullptr;   // только SIM-эксперименты: время смерти по PID (команда), -1 — не умрёт
    int curCmd = 0;

    double kappa(int f) const {
        if (coTot[f] < P[21]) return 0;
        double b0 = birthV[f][0], b1 = birthV[f][1], bt = b0 + b1;
        if (bt <= 0) return 0;
        if (std::min(b0, b1) < P[36] * bt) return 0;   // признак почти константа — каппа вырождена
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
        // гистерезис: текущий признак меняем, только если другой лучше с запасом
        int bf = -1; double bs = -1, cs = -1;
        for (int f = 0; f < NF; f++) {
            double k = kappa(f);
            if (k < P[22]) continue;
            double s = k + std::fabs(pdie(f, 0) - pdie(f, 1));
            if (f == curFeat) cs = s;
            if (s > bs) { bs = s; bf = f; }
        }
        if (bf < 0) { curFeat = -1; return -1; }
        if (cs >= 0 && bs < cs + 0.2) bf = curFeat;
        curFeat = bf;
        int v = (fv >> bf) & 1;
        // тот, кто умирает чаще, — в поток B (правый край, у границы роста)
        double pv = pdie(bf, v), po = pdie(bf, v ^ 1);
        if (std::fabs(pv - po) < 0.1) return v;
        return pv > po ? 1 : 0;
    }

    bool kin(int p, int st) const { return p > 0 && stream[p] == st && runOf[p] == curRun; }

    // положить в промежуток gs длины s: выбрать сторону, вернуть D
    int chooseSide(int gs, int s, int L, int st) {
        int r = s - L;
        if (r == 0) return gs;
        int lp = pidLeftOf(gs), rp = pidRightOf(gs + s);
        if (st >= 0) {
            if (P[29] != 0) {
                bool kl = kin(lp, st), kr = kin(rp, st);
                if (kl && !kr) return gs;
                if (kr && !kl) return gs + r;
            }
            // прижимаем к соседу из того же потока
            int lc = lp ? stream[lp] : -2, rc = rp ? stream[rp] : -2;
            if (lc == st && rc != st) return gs;
            if (rc == st && lc != st) return gs + r;
            return st == 1 ? gs + r : gs;
        }
        if (P[8] == 4) { if (hz2(lp) > hz2(rp)) return gs + r; return gs; }
        if (P[8] == 1) { if (lp > rp) return gs + r; }
        else if (P[8] == 3) return gs + r;
        return gs;
    }

    std::pair<ll,ll> place(int L) {
        int m = minRecent();
        int pid = (int)pos.size();
        if (L == lastL) runIdx++; else { runIdx = 0; curRun++; }
        runOf.push_back(curRun);
        lastL = L;
        feat.push_back((unsigned char)((runIdx & 1) | ((pid & 1) << 1)));
        int st = chooseStream(pid, L);
        stream.push_back((signed char)st);
        bVol = bVol * P[27] + (st == 1 ? L : 0);
        if (st >= 0) bornLab[st] += 1;
#ifdef SIM
        if (getenv("SFDBG") && pid % atoi(getenv("SFDBG")) == 0)
            fprintf(stderr, "pid=%d L=%d feat=%d st=%d k0=%.3f k1=%.3f tot=%.1f,%.1f F=%d T=%lld\n", pid, L, feat[pid], st, kappa(0), kappa(1), coTot[0], coTot[1], F, T);
#endif
        for (int f = 0; f < NF; f++) {
            birthV[f][0] *= P[23]; birthV[f][1] *= P[23];
            deathV[f][0] *= P[23]; deathV[f][1] *= P[23];
            birthV[f][(feat[pid] >> f) & 1] += 1;
        }
        if (F > T) F = (int)T;
        int ts = tailStart();
        if (P[14] == 5 && ++wAge >= 64) rebuildW();
        double selfHz = 0;
        if (P[18] == 4 && recDelTot >= 1)
            selfHz = std::min(1.0, (recDel[L] / recDelTot) / ((double)std::max(1, aliveCntL[L]) / std::max(1, aliveTot)) / P[30]);
        double bestC = 1e18; int bgs = -1, bs = 0;
        int M = P[4] > 0 ? (int)P[4] : 1;  // сколько промежутков каждой длины смотреть
        auto eval = [&](int gs, int s) {
            double cst;
            if (P[13] > 0 && T < P[13]) {
                // крошечная траншея: точное попадание (если сосед не «умирающий»), иначе самая большая дыра
                if (s == L) cst = -1000 + P[32] * std::max(pdieLen(pidLeftOf(gs)), pdieLen(pidRightOf(gs + s)));
                else cst = -s;
            } else if (P[14] == 1) {
                cst = gs;                                   // first fit по адресу
            } else if (P[14] == 2) {
                cst = -s;                                   // worst fit
            } else if (P[14] == 3) {
                cst = s == L ? -1e9 : -s;                   // точное, иначе worst
            } else if (P[14] == 5) {
                // изменение ожидаемых потерь + разрыв ничьих по best fit
                cst = Wf(s - L) - Wf(s) + 1e-3 * (s - L);
            } else if (P[14] == 4) {
                cst = -(double)gs;                          // last fit
            } else {
                cst = s == L ? -P[5] : fragCost(s - L, m) + P[6] * (s - L);
                if (P[17] != 0 && s == L) {
                    // точное попадание «замуровывает» дыру у соседа, который скоро умрёт (дыра бы выросла)
                    double h;
                    if (P[39] > 0) {
                        // дыра вырастет после смерти соседа до s + len(соседа) — полезно, только если такая длина
                        // часто встречается среди недавних укладок
                        h = 0;
                        int nb[2] = {pidLeftOf(gs), pidRightOf(gs + s)};
                        for (int k = 0; k < 2; k++) {
                            if (nb[k] <= 0) continue;
                            int mg = s + len[nb[k]];
                            if (mg > MAXL || cnt[mg] < P[39] * cntTot) continue;
                            h = std::max(h, hz2(nb[k]));
                        }
                    } else h = std::max(hz2(pidLeftOf(gs)), hz2(pidRightOf(gs + s)));
                    if (P[18] == 4) h -= selfHz;   // вредно, только если новая труба живёт дольше соседа
                    if (h >= P[28]) cst += P[17] * h;
                }
                if (P[16] != 0 && s > L && s - L < m) {
                    // мёртвый остаток оживёт, только если умрёт сосед: штраф ~ r * P(оба соседа выживут)
                    double hl = hazard(pidLeftOf(gs)), hr = hazard(pidRightOf(gs + s));
                    cst += P[16] * (s - L) * (1 - std::max(hl, hr));
                }
            }
            if (P[29] != 0 && st >= 0 && (kin(pidLeftOf(gs), st) || kin(pidRightOf(gs + s), st))) cst -= P[29];
            if (P[40 - 1 + 0] != 0 && false) {}
            if (P[15] > 0 && P[33] == 2 && st < 0 && s > L && s - L < m && s - L >= P[15]) cst = 1e17;   // мёртвый остаток — не кандидат
            if (tieAge != 0) {
                // разбивка ничьих по возрасту соседей (молодые = большой PID)
                int lp = pidLeftOf(gs), rp = pidRightOf(gs + s);
                double young = std::max(lp, rp), old = std::min(lp, rp);
                if (s == L) cst += tieAge * 1e-7 * old;        // точное — между старыми
                else cst -= tieAge * 1e-7 * young;             // остаток — к молодому
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
#ifdef SIM
        if (forcePid == pid) {
            // принудительное отклонение: forceAlt = -1 — хвост, иначе дыра размера forceAlt
            if (forceAlt < 0) bgs = -1;
            else { bgs = -1; for (int g : bucket[forceAlt]) if (g != ts && g + forceAlt != F) { bgs = g; bs = forceAlt; bestC = 0; break; } }
        }
        if (listPid == pid) {
            optBest = bgs >= 0 ? bs : -1; optList.clear();
            for (int s = L; s < NB; s++) for (int g : bucket[s]) { if (g == ts || g + s == F) continue; optList.push_back(s); break; }
            optList.push_back(-1);
        }
        if (epsRand > 0 && (rng() % 1000) < epsRand) {
            // случайный выбор: любая подходящая дыра (по одной на размер) или хвост
            int opts[NB + 2], no = 0;
            for (int s = L; s < NB; s++) for (int g : bucket[s]) { if (g == ts || g + s == F) continue; opts[no++] = s; break; }
            int pick = (int)(rng() % (no + 1));
            if (getenv("SFEPSLOG")) fprintf(stderr, "dev t=%d pid=%d L=%d best=%s%d pick=%s%d opts=%d\n", curCmd, pid, L, bgs >= 0 ? "hole" : "tail", bgs >= 0 ? bs : 0, pick == no ? "tail" : "hole", pick == no ? 0 : opts[pick], no);
            if (pick == no) bgs = -1;
            else { int s = opts[pick]; for (int g : bucket[s]) if (g != ts && g + s != F) { bgs = g; bs = s; break; } bestC = 0; }
        }
#endif
        ll N = 0; int D;
        bool useGap = bgs >= 0 && bestC < P[7] + 1e17 * (P[7] == 0);
        if (st == 1 && P[12] == 2) useGap = false;
        if (st == 1 && P[12] == 3 && bornLab[1] >= 20 && bornLab[0] >= 20 &&
            diedLab[1] / bornLab[1] - diedLab[0] / bornLab[0] > P[9]) useGap = false;   // прогноз подтвердился: B — к границе роста
        if (useGap && P[14] == 5 && bestC > P[19]) useGap = false;   // дыра не стоит того — в хвост
        // P15: если лучший промежуток оставляет мёртвый остаток >= P15 — лучше в хвост
        if (bgs >= 0 && bestC >= 1e17) useGap = false;
        if (useGap && P[15] > 0 && P[33] != 2 && !(P[33] == 1 && st >= 0)) {
            int r = bs - L;
            if (r > 0 && r < m) {
                double keep = P[25] ? 1 - std::max(hz2(pidLeftOf(bgs)), hz2(pidRightOf(bgs + bs))) : 1;
                if (r * keep >= (P[34] > 0 ? P[34] * m : P[15])) useGap = false;
            }
        }
        if (useGap) {
            D = chooseSide(bgs, bs, L, st);
            takeFrom(bgs, D, L);
        } else {
            int fgs = (F > 0 && tag[F] < 0) ? -tag[F] - 1 : F;
            int fl = F - fgs;
            int chunk = (int)P[20];
            if (P[26] > 0) chunk = std::min(chunk, std::max(L, (int)(P[26] * bVol)));
            if (P[35] > 0) chunk = std::min(chunk, std::max(2 * L, (int)(P[35] * T)));   // не больше доли траншеи
            if (P[37] > 0 && st == 1) chunk = std::min(chunk, std::max(L, (int)(P[37] * fbVol)));
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
#ifdef SIM
        if (getenv("SFPL") && pid >= atoi(getenv("SFPL")) && pid < atoi(getenv("SFPL")) + 30)
            fprintf(stderr, "place pid=%d L=%d st=%d D=%d N=%lld gap=%s bs=%d F=%d T=%lld\n", pid, L, st, D, N, useGap ? "hole" : "front", bs, F, T);
#endif
        for (int i = 0; i < L; i++) owner[D + i] = pid;
        fbVol = fbVol * P[38] + ((st == 1 && !useGap) ? L : 0);
        pos.push_back(D); len.push_back(L); alive.push_back(1); birthL[L] += 1; aliveCntL[L]++; aliveTot++;
        if ((int)aliveBit.size() <= pid) aliveBit.resize(std::max<size_t>(pid + 1, aliveBit.size() * 2 + 16), 0), rebuildBit();
        bitAdd(pid, 1);
        observe(L); nPlaced++;
        lastWasDel = false;
        return {N, D};
    }
    void rebuildBit() {
        // после расширения массива Фенвика пересобрать по живым
        std::fill(aliveBit.begin(), aliveBit.end(), 0);
        for (int p = 1; p < (int)pos.size() && p < (int)aliveBit.size(); p++) if (alive[p]) bitAdd(p, 1);
    }
    std::vector<char> alive;
    ll remove(int pid) {
        int d = pos[pid], l = len[pid];
        alive[pid] = 0; bitAdd(pid, -1); deathL[l] += 1;
        for (int k = 1; k <= MAXL; k++) recDel[k] *= P[31];
        recDelTot = recDelTot * P[31] + 1; recDel[l] += 1; aliveCntL[l]--; aliveTot--;
        deathAge[ageBucket((int)pos.size() - 1 - pid)] += 1;
        int e = freeRange(d, l);
        if (d == F) F = e;
        unsigned char fv = feat[pid];
        if (stream[pid] >= 0) diedLab[stream[pid]] += 1;
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
    Alloc() { pos.push_back(-1); len.push_back(0); feat.push_back(0); stream.push_back(-1); alive.push_back(0); runOf.push_back(-1); ensure(1024); }
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
// T в конце прогона (или -1 при ошибке); force/list — для анализа важности решений
static ll runOnce(const std::vector<int>& cmd, int forcePid, int forceAlt, int listPid, Alloc** keep) {
    Alloc* A = new Alloc();
    A->forcePid = forcePid; A->forceAlt = forceAlt; A->listPid = listPid;
    int ci = 0;
    for (int c : cmd) { A->curCmd = ci++; if (c > 0) A->place(c); else A->remove(-c); }
    ll T = A->T;
    if (keep) *keep = A; else delete A;
    return T;
}
int main(int argc, char** argv) {
    if (const char* s = getenv("SF")) {
        // формат: "i=v i=v ..."
        char* p = (char*)s;
        while (*p) {
            char* q; long i = strtol(p, &q, 10); if (q == p || *q != '=') break;
            p = q + 1; double v = strtod(p, &q); if (q == p) break; if (i >= 0 && i < 40) P[i] = v; p = q;
        }
    }
    bool quiet = getenv("SFQ") != nullptr;
    if (getenv("SFIMPORT")) {
        double w; std::vector<int> cmd = readTest(argv[1], w);
        ll T0 = runOnce(cmd, -1, 0, -1, nullptr);
        int np = 0; for (int c : cmd) if (c > 0) np++;
        printf("base T=%lld\n", T0);
        for (int p = 1; p <= np; p++) {
            Alloc* A; runOnce(cmd, -1, 0, p, &A);
            std::vector<int> opts = A->optList; int best = A->optBest; int L = A->len[p];
            delete A;
            for (int o : opts) {
                if (o == best) continue;
                ll T = runOnce(cmd, p, o, -1, nullptr);
                if (T != T0) printf("pid=%d L=%d best=%d alt=%d dT=%lld\n", p, L, best, o, T - T0);
            }
        }
        return 0;
    }
    double total = 0, wsum = 0;
    for (int a = 1; a < argc; a++) {
        double w; std::vector<int> cmd = readTest(argv[a], w);
        Alloc A;
        std::vector<int> death(1, -1);
        { int p = 0; for (size_t i = 0; i < cmd.size(); i++) { if (cmd[i] > 0) { death.push_back(-1); p++; } else death[-cmd[i]] = (int)i; } }
        A.oracleDeath = death.data();
        if (getenv("SFTIE")) A.tieAge = atof(getenv("SFTIE"));
        if (getenv("SFFORCE")) { A.forcePid = atoi(getenv("SFFORCE")); A.forceAlt = atoi(getenv("SFALT")); }
        if (getenv("SFLIST")) A.listPid = atoi(getenv("SFLIST"));
        if (getenv("SFEPS")) { A.epsRand = atoi(getenv("SFEPS")); A.rs ^= (unsigned long long)atoll(getenv("SFSEED")) * 0x9E3779B97F4A7C15ULL + 1; }
        std::vector<int> cells;  // pid или 0
        std::vector<int> lastLen, freedAt;
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
                cells.resize(T, 0); lastLen.resize(T, 0); freedAt.resize(T, -1);
                for (int i = 0; i < c; i++) if (cells[r.second + i]) { ok = false; break; }
                if (!ok) break;
                int pid = (int)A.pos.size() - 1;
                if (getenv("SFDIG") && r.first > 0) {
                    ll holes, small; A.holeStats(holes, small);
                    int mx = 0; for (int s = NB - 1; s > 0; s--) if (!A.bucket[s].empty()) { for (int g : A.bucket[s]) if (g != A.tailStart() && g + s != A.F) { mx = s; break; } if (mx) break; }
                    printf("  dig t=%lld L=%d N=%lld T=%lld occ=%lld holes=%lld maxhole=%d\n", ci, c, r.first, A.T, occ, holes, mx);
                }
                for (int i = 0; i < c; i++) cells[r.second + i] = pid, lastLen[r.second + i] = c;
                occ += c;
            } else {
                int pid = -c;
                ll d = A.remove(pid);
                if (d < 0 || d >= (ll)cells.size() || cells[d] != pid) { ok = false; break; }
                int l = A.len[pid];
                for (int i = 0; i < l; i++) cells[d + i] = 0, freedAt[d + i] = (int)ci;
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
        if (getenv("SFFREED") && ok) {
            int from = atoi(getenv("SFFREED")); long long n = (long long)cmd.size();
            long long b[8] = {0}; // освобождено: <100, <300, <1000, <3000, <10000, <30000, раньше, никогда
            int lim[6] = {100, 300, 1000, 3000, 10000, 30000};
            for (size_t i = from; i < cells.size(); i++) if (!cells[i]) {
                if (freedAt[i] < 0) { b[7]++; continue; }
                long long age = n - freedAt[i]; int k = 0; while (k < 6 && age >= lim[k]) k++; b[k]++;
            }
            printf("  freed-age from %d: <100:%lld <300:%lld <1k:%lld <3k:%lld <10k:%lld <30k:%lld older:%lld never:%lld\n", from, b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]);
        }
        if (getenv("SFMAP") && ok) {
            // карта траншеи: длины труб (с возрастом) и дыры [h]
            size_t i = 0; int np = (int)A.pos.size() - 1;
            while (i < cells.size()) {
                if (cells[i]) { int p = cells[i]; printf("%d", A.len[p]); if (getenv("SFMAPAGE")) printf("(%d)", np - p); printf(" "); i += A.len[p]; }
                else { size_t j = i; while (j < cells.size() && !cells[j]) j++; printf("[%zu] ", j - i); i = j; }
            }
            printf("\n");
        }
        if (getenv("SFHOLES") && ok) {
            int from = atoi(getenv("SFHOLES"));
            long long h[200] = {0}; long long tot = 0;
            for (size_t i = from; i < cells.size();) { if (cells[i]) { i++; continue; } size_t j = i; while (j < cells.size() && !cells[j]) j++; int s = (int)std::min<size_t>(j - i, 199); if (j < cells.size()) { h[s] += j - i; tot += j - i; } i = j; }
            printf("  holes from %d: total %lld; by size (meters):", from, tot);
            for (int s = 1; s < 200; s++) if (h[s]) printf(" %d:%lld", s, h[s]);
            printf("\n");
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
