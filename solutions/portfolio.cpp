// Вариант B: портфель политик размещения + выбор политики по истории (теневые симуляции).
// Идея: набор «чистых» политик (best/first/last/worst/exact fit, арены для скоротечных труб по предсказателю времени жизни)
// прогоняется параллельно в теневых копиях на всей истории команд; реальная траншея ведётся текущим лидером.
// Сборка решения:  c++ -std=c++17 -O2 -o build/portfolio solutions/portfolio.cpp
// Стенд (внутрипроцессный прогон тестов): scripts/sim.cpp, сборка scripts/build_sim.bat; эксперименты: scripts/exp_*.cpp
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
#ifdef _MSC_VER
#include <intrin.h>
#endif
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#else
#include <time.h>
#endif

// CPU-время процесса в секундах (для аварийного упрощения при медленной проверяющей машине)
static double cpuSeconds() {
#ifdef _WIN32
    FILETIME c, e, k, u;
    GetProcessTimes(GetCurrentProcess(), &c, &e, &k, &u);
    auto t = [](const FILETIME& f) { return (double)(((uint64_t)f.dwHighDateTime << 32) | f.dwLowDateTime) * 1e-7; };
    return t(k) + t(u);
#else
    timespec ts;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
#endif
}

static inline int ctz64(uint64_t x) {
#ifdef _MSC_VER
    unsigned long r;
    _BitScanForward64(&r, x);
    return (int)r;
#else
    return __builtin_ctzll(x);
#endif
}
static inline int clz64idx(uint64_t x) {  // индекс старшего бита
#ifdef _MSC_VER
    unsigned long r;
    _BitScanReverse64(&r, x);
    return (int)r;
#else
    return 63 - __builtin_clzll(x);
#endif
}

#ifndef LEAN_PER_CMD
#define LEAN_PER_CMD 21e-6   // порог CPU на команду для аварийного режима: 5 с / 200000 = 25 мкс, берём запас
#endif

// ---------------------------------------------------------------- политики

enum Sel { BEST, FIRST, LAST_UNUSED, WORST, EXACTWORST, TFIT, SCORED, WFN, XFIRST };   // LAST убран (кучи хранят только минимальный адрес)

struct Pol {
    const char* name;
    Sel sel;          // как выбирать дыру для обычных (P) труб
    bool tieHigh;     // (не используется; оставлено для совместимости пулов)
    bool right;       // класть к правому краю дыры
    bool tailFirst;   // хвост — обычный кандидат: если влезает, класть туда до поиска дыр
    int js;           // горизонт предсказателя для класса T (скоротечные); -1 — без арен, все трубы класса P
    int Cdiv, Cmin, Cmax;  // размер чанка арены: clamp(trench / Cdiv, Cmin, Cmax)
    int Tsmall;       // P-трубы могут занимать T-дыры размером <= Tsmall
    int tuU = 0;      // 1: T-трубы тоже могут занимать U-дыры (best fit по обеим)
    int Cmode = 0;    // 0: чанк = Cdiv * (средняя длина T-труб); 1: чанк = liveT / Cdiv
    int th = 0;       // порог остатка для TFIT
    float tau = 0.5f; // класс T, если предсказанная вероятность смерти >= tau
    float sw0 = 0, swsmall = 0, swsize = 0;  // SCORED/WFN: бонус точного попадания, штраф за малый остаток (< th), вес размера дыры
    float wkap = 1e9f;                       // WFN: если лучшая стоимость > wkap * L — труба идёт в хвост
};

static const Pol POLS_INIT[] = {
    //  имя                        sel         tieHi  right  tailF  js  Cdiv Cmin Cmax Tsm tuU Cmode th  tau
    {"bf-L",                      BEST,       false, false, false, -1,  0,   0,   0,   0,  0,  0,    0, 0.5f},   // 0: основная стартовая (безопасная)
    {"a99_E_d25_c100_t0_u1",      EXACTWORST, false, false, false, 99, 25, 100, 4000,  0,  1,  0,    0, 0.7f},   // арены, горизонт подбирается адаптивно
    {"ff-L",                      FIRST,      false, false, false, -1,  0,   0,   0,   0,  0,  0,    0, 0.5f},
    {"xw-R",                      EXACTWORST, false, true,  false, -1,  0,   0,   0,   0,  0,  0,    0, 0.5f},   // первые startCmds команд и точные попадания
    {"a99_B_d12_c40_t40_u1",      BEST,       false, false, false, 99, 12,  40, 4000, 40,  1,  0,    0, 0.7f},
    {"a99_B_d6_c40_t40_u0",       BEST,       false, false, false, 99,  6,  40, 4000, 40,  0,  0,    0, 0.7f},
};
static const int NPOL0 = sizeof(POLS_INIT) / sizeof(POLS_INIT[0]);
// Таблица политик изменяема, чтобы стенд мог перебирать параметры (в бою не меняется).
static std::vector<Pol>& polTable() {
    static std::vector<Pol> v(POLS_INIT, POLS_INIT + NPOL0);
    return v;
}
#define POLS (polTable())
#define NPOL ((int)polTable().size())

// ---------------------------------------------------------------- хеш-таблица дыр и множество свободных дыр одного тега

// Открытая адресация, ключ — адрес (>= 0), значение — длина и тег; удаление со сдвигом назад (без надгробий).
struct HTab {
    std::vector<int> key, val;
    std::vector<uint8_t> tag;
    unsigned cap = 0, msk = 0, count = 0;
    HTab() { init(256); }
    void init(unsigned c) { cap = c; msk = c - 1; key.assign(c, -1); val.assign(c, 0); tag.assign(c, 0); count = 0; }
    static inline unsigned hsh(int k) { return (unsigned)k * 2654435761u; }
    inline unsigned home(int k) const { return (hsh(k) >> 7) & msk; }
    inline int find(int k) const {   // индекс слота или -1
        unsigned i = home(k);
        while (key[i] != -1) { if (key[i] == k) return (int)i; i = (i + 1) & msk; }
        return -1;
    }
    void put(int k, int v, int t) {
        if ((count + 1) * 2 > cap) grow();
        unsigned i = home(k);
        while (key[i] != -1) i = (i + 1) & msk;
        key[i] = k; val[i] = v; tag[i] = (uint8_t)t; count++;
    }
    void grow() {
        std::vector<int> ok = std::move(key), ov = std::move(val);
        std::vector<uint8_t> ot = std::move(tag);
        unsigned oc = cap;
        init(cap * 2);
        for (unsigned i = 0; i < oc; i++) if (ok[i] != -1) put(ok[i], ov[i], ot[i]);
    }
    void erase(int k) {
        int ii = find(k);
        unsigned i = (unsigned)ii, j = i;
        for (;;) {
            j = (j + 1) & msk;
            if (key[j] == -1) break;
            unsigned h = home(key[j]);
            if (i <= j ? (i < h && h <= j) : (i < h || h <= j)) continue;
            key[i] = key[j]; val[i] = val[j]; tag[i] = tag[j];
            i = j;
        }
        key[i] = -1;
        count--;
    }
};

// Дыры одного тега: корзины по размеру 1..39 — минимальные кучи адресов с ленивым удалением (актуальность проверяется
// по таблице дыр), большие дыры (>= 40) — упорядоченные множества. Счётчики cnt[] точные.
struct HoleSet {
    std::vector<int> heap[40];            // min-кучи адресов дыр длины s (могут содержать устаревшие записи)
    int cnt[40] = {0};                    // число актуальных дыр длины s (1..39)
    std::set<std::pair<int, int>> big;    // дыры длины >= 40: (длина, адрес)
    std::set<int> bigAddr;                // они же по адресу
    uint64_t mask = 0;                    // бит s (1..39) — есть дыра длины s; бит 40 — есть большая
    long long bigSum = 0;                 // суммарная длина больших дыр

    void add(int a, int h) {
        if (h < 40) {
            heap[h].push_back(a);
            std::push_heap(heap[h].begin(), heap[h].end(), std::greater<int>());
            cnt[h]++; mask |= 1ull << h;
        } else { big.insert({h, a}); bigAddr.insert(a); mask |= 1ull << 40; bigSum += h; }
    }
    void del(int a, int h) {
        if (h < 40) { cnt[h]--; if (!cnt[h]) mask &= ~(1ull << h); }
        else { big.erase({h, a}); bigAddr.erase(a); bigSum -= h; if (big.empty()) mask &= ~(1ull << 40); }
    }
    // суммарная длина дыр короче q (непригодных для труб длины >= q)
    long long unusable(int q) const {
        long long u = 0;
        for (int sI = 1; sI < q && sI < 40; sI++) u += (long long)sI * cnt[sI];
        return u;
    }
};

// ---------------------------------------------------------------- политики: выбор дыры

// наименьший адрес актуальной дыры длины s в наборе hs (tg — тег набора); очищает устаревшие верхушки кучи
static const int WMAX = 80;
static inline int topAddr(HoleSet& hs, const HTab& st, int s, int tg) {
    std::vector<int>& hp = hs.heap[s];
    for (;;) {
        int a = hp.front();
        int i = st.find(a);
        if (i >= 0 && st.val[i] == s && st.tag[i] == tg) return a;
        std::pop_heap(hp.begin(), hp.end(), std::greater<int>());
        hp.pop_back();
    }
}

// ---------------------------------------------------------------- состояние траншеи

struct Alloc {
    int trench = 0;      // длина траншеи
    int tailStart = 0;   // начало свободного хвоста (хвост = [tailStart, trench))
    HTab startTab;       // дыры по началу: начало -> (длина, тег); дыры внутри траншеи, не касаются конца
    HTab endTab;         // дыры по концу: конец -> длина
    HoleSet hs[2];       // 0 — обычные (U), 1 — арены/дыры скоротечных (T)
    std::vector<int> pos, len;            // по PID (с 1)
    std::vector<uint8_t> cls;             // класс трубы: 0 = P (обычная), 1 = T (скоротечная)
    long long live = 0;
    long long liveT = 0;                  // суммарная длина живых труб класса T
    double emaT = 0;                      // скользящее среднее длины недавно уложенных T-труб
    const Pol* curPol = nullptr;
    const double* W = nullptr;            // таблица «ожидаемого непригодного остатка» W[s], s = 0..WMAX (для WFN)

    Alloc() { pos.assign(1, -1); len.assign(1, 0); cls.assign(1, 0); }

    void addHole(int a, int h, int tag) { startTab.put(a, h, tag); endTab.put(a + h, h, 0); hs[tag].add(a, h); }
    void delHole(int a) {
        int i = startTab.find(a);
        int h = startTab.val[i], tg = startTab.tag[i];
        hs[tg].del(a, h);
        endTab.erase(a + h);
        startTab.erase(a);
    }
    int holeLen(int a) const { return startTab.val[startTab.find(a)]; }
    template <class F> void forEachHole(F f) const {   // f(начало, длина, тег)
        for (unsigned i = 0; i < startTab.cap; i++) if (startTab.key[i] != -1) f(startTab.key[i], startTab.val[i], (int)startTab.tag[i]);
    }

    // вырезать [D, D+L) из дыры (a, h, tag)
    void carve(int a, int h, int tag, int D, int L) {
        delHole(a);
        if (D > a) addHole(a, D - a, tag);
        if (D + L < a + h) addHole(D + L, a + h - (D + L), tag);
    }
    void commit(int D, int L, int c) {
        pos.push_back(D); len.push_back(L); cls.push_back((uint8_t)c); live += L;
        if (c) { liveT += L; emaT = emaT == 0 ? L : emaT * 0.9 + L * 0.1; }
    }

    // кандидат из набора tg по правилу; maxSize ограничивает размер дыры. false — нет подходящей. h получает длину дыры.
    bool pick(int tg, Sel sel, int L, int maxSize, int th, int& a, int& h) {
        HoleSet& H = hs[tg];
        uint64_t mm = H.mask & (~0ull << L);
        if (maxSize < 40) mm &= (2ull << maxSize) - 1;
        else if ((mm >> 40 & 1) && maxSize < (1 << 30) && H.big.begin()->first > maxSize) mm &= ~(1ull << 40);
        if (!mm) return false;
        switch (sel) {
        case BEST: {
            int s = ctz64(mm);
            if (s < 40) { a = topAddr(H, startTab, s, tg); h = s; }
            else { a = H.big.begin()->second; h = H.big.begin()->first; }
            return true;
        }
        case XFIRST:
            if (L < 40 && (H.mask >> L & 1) && L <= maxSize) { a = topAddr(H, startTab, L, tg); h = L; return true; }
            /* fallthrough */
        case FIRST: {
            int best = 0x7fffffff, bh = 0;
            for (uint64_t m = mm & ((1ull << 40) - 1); m; m &= m - 1) {
                int s = ctz64(m);
                int x = topAddr(H, startTab, s, tg);
                if (x < best) { best = x; bh = s; }
            }
            if (mm >> 40 & 1) {
                int x = *H.bigAddr.begin();
                if (x < best) { best = x; bh = holeLen(x); }
            }
            a = best; h = bh;
            return true;
        }
        case WORST:
            if (mm >> 40 & 1) { a = H.big.rbegin()->second; h = H.big.rbegin()->first; }
            else { int s = clz64idx(mm); a = topAddr(H, startTab, s, tg); h = s; }
            return true;
        case EXACTWORST:
            if (L < 40 && (H.mask >> L & 1) && L <= maxSize) { a = topAddr(H, startTab, L, tg); h = L; return true; }
            if (mm >> 40 & 1) { a = H.big.rbegin()->second; h = H.big.rbegin()->first; }
            else { int s = clz64idx(mm); a = topAddr(H, startTab, s, tg); h = s; }
            return true;
        case SCORED: {   // минимум стоимости: -sw0 при точном попадании + swsmall*(th-r) при 0<r<th + swsize*размер
            const Pol* pp = curPol;
            double bestC = 1e300; int bs = -1;
            for (uint64_t m = mm & ((1ull << 40) - 1); m; m &= m - 1) {
                int sI = ctz64(m), r = sI - L;
                double c = pp->swsize * sI;
                if (r == 0) c -= pp->sw0;
                else if (r < th) c += pp->swsmall * (th - r);
                if (c < bestC) { bestC = c; bs = sI; }
            }
            if (mm >> 40 & 1) {
                for (int k = 0; k < 2; k++) {
                    int sI = k == 0 ? H.big.begin()->first : H.big.rbegin()->first;
                    double c = pp->swsize * sI;   // остаток большой
                    if (c < bestC) { bestC = c; bs = sI; }
                }
            }
            if (bs < 40) { a = topAddr(H, startTab, bs, tg); h = bs; }
            else { auto it = H.big.lower_bound({bs, -1}); a = it->second; h = it->first; }
            return true;
        }
        case WFN: {   // минимум прироста ожидаемых потерь: W[h-L] - W[h] (+ предпочтение размера, бонус точного попадания)
            const Pol* pp = curPol;
            double bestC = 1e300; int bs = -1;
            for (uint64_t m = mm & ((1ull << 40) - 1); m; m &= m - 1) {
                int sI = ctz64(m);
                double c = W[sI - L] - W[sI] + pp->swsize * sI;
                if (sI == L) c -= pp->sw0;
                if (c < bestC) { bestC = c; bs = sI; }
            }
            if (mm >> 40 & 1) {
                for (int k = 0; k < 2; k++) {
                    int sI = k == 0 ? H.big.begin()->first : H.big.rbegin()->first;
                    int s2 = std::min(sI, WMAX), r2 = std::min(sI - L, WMAX);
                    double c = W[r2] - W[s2] + pp->swsize * sI;
                    if (c < bestC) { bestC = c; bs = sI; }
                }
            }
            if (bestC > pp->wkap * L) return false;
            if (bs < 40) { a = topAddr(H, startTab, bs, tg); h = bs; }
            else { auto it = H.big.lower_bound({bs, -1}); a = it->second; h = it->first; }
            return true;
        }
        case TFIT: {
            if (L < 40 && (H.mask >> L & 1) && L <= maxSize) { a = topAddr(H, startTab, L, tg); h = L; return true; }
            uint64_t hi = (L + th < 64) ? (mm & (~0ull << (L + th))) : 0;
            int s = hi ? ctz64(hi) : clz64idx(mm);
            if (s < 40) { a = topAddr(H, startTab, s, tg); h = s; }
            else { a = H.big.begin()->second; h = H.big.begin()->first; }
            return true;
        }
        default: break;
        }
        return false;
    }

    // уложить трубу длины L класса c; возвращает N (докопка в 20-метровых сегментах) и D (позиция)
    void place(const Pol& p, int L, int c, int& N, int& D) {
        N = 0;
        curPol = &p;
        int a, h;
        if (c == 1) {  // скоротечная: best fit по T-дырам (и, если разрешено, по U-дырам), иначе новый чанк в хвосте
            int a1 = 0, h1 = 0, a0 = 0, h0 = 0;
            bool f1 = pick(1, BEST, L, 1 << 30, 0, a1, h1);
            bool f0 = p.tuU && pick(0, BEST, L, 1 << 30, 0, a0, h0);
            if (f1 || f0) {
                bool useT = f1 && (!f0 || h1 <= h0);
                a = useT ? a1 : a0; h = useT ? h1 : h0;
                D = a; carve(a, h, useT ? 1 : 0, D, L); commit(D, L, 1); return;
            }
            int C = std::min(std::max((int)(p.Cmode == 0 ? emaT * p.Cdiv : liveT / p.Cdiv), p.Cmin), p.Cmax);
            int need = std::max(C, L);
            int tailLen = trench - tailStart;
            if (tailLen < need) { N = (need - tailLen + 19) / 20; trench += 20 * N; }
            D = tailStart;
            tailStart += need;
            if (need > L) addHole(D + L, need - L, 1);
            commit(D, L, 1);
            return;
        }
        int tailLen = trench - tailStart;
        if (p.tailFirst && tailLen >= L) {
            D = tailStart; tailStart += L; commit(D, L, 0); return;
        }
        // P: кандидаты из U-дыр и малых T-дыр
        int a0 = 0, h0 = 0, a1 = 0, h1 = 0;
        bool f0 = pick(0, p.sel, L, 1 << 30, p.th, a0, h0);
        bool f1 = pick(1, p.sel, L, p.js >= 0 ? p.Tsmall : (1 << 30), p.th, a1, h1);   // политики без арен видят T-дыры как обычные
        if (f0 || f1) {
            bool useT;
            if (!f1) useT = false;
            else if (!f0) useT = true;
            else {
                switch (p.sel) {
                case BEST: useT = h1 < h0 || (h1 == h0 && a1 < a0); break;
                case FIRST: useT = a1 < a0; break;
                default: useT = h1 > h0; break;
                }
            }
            int tag = useT ? 1 : 0;
            a = useT ? a1 : a0; h = useT ? h1 : h0;
            if (p.right) D = a + h - L; else D = a;
            carve(a, h, tag, D, L);
            commit(D, L, 0);
            return;
        }
        if (tailLen < L) { N = (L - tailLen + 19) / 20; trench += 20 * N; }
        D = tailStart;
        tailStart += L;
        commit(D, L, 0);
    }

    int remove(int pid) {
        int p = pos[pid], l = len[pid];
        int tag = cls[pid];
        pos[pid] = -1;
        live -= l;
        if (cls[pid]) liveT -= l;
        int s = p, e = p + l;
        int li = endTab.find(p);   // дыра, заканчивающаяся в p
        if (li >= 0) {
            int hl = endTab.val[li];
            int st = p - hl;
            tag = std::max(tag, (int)startTab.tag[startTab.find(st)]);
            s = st;
            delHole(st);
        }
        if (e == tailStart) {
            tailStart = s;
        } else {
            int ri = startTab.find(e);   // дыра, начинающаяся в e
            if (ri >= 0) { tag = std::max(tag, (int)startTab.tag[ri]); int hl = startTab.val[ri]; delHole(e); e += hl; }
            addHole(s, e - s, tag);
        }
        return p;
    }
};

// ---------------------------------------------------------------- предсказатель времени жизни
// Для каждой трубы при укладке строятся ключи контекстов (семейства признаков). По завершённым наблюдениям
// («умерла за H команд» / «дожила до возраста H») каждое семейство ведёт таблицу P(смерть за H | ключ) по набору горизонтов.
// Качество семейства оценивается онлайн (префиксная кросс-энтропия с затуханием); для прогноза берётся лучшее.
// Пространства ключей маленькие и плотные — таблицы без хеширования.
static inline float fastLog(float x) {   // натуральный логарифм, погрешность ~1e-3 (достаточно для сравнения семейств)
    uint32_t i;
    std::memcpy(&i, &x, sizeof i);
    float e = (float)((int)(i >> 23) - 127);
    i = (i & 0x007FFFFF) | 0x3F800000;       // мантисса в [1,2)
    float m;
    std::memcpy(&m, &i, sizeof m);
    float lg2 = e + (-0.34484843f * m + 2.02466578f) * m - 0.67487759f;
    return lg2 * 0.69314718f;
}

struct Predictor {
    static constexpr int NH = 8;              // горизонты H_j = 8 * 4^j команд (до 131072)
    static constexpr int NF = 15;             // семейства признаков
    static int H(int j) { return 8 << (2 * j); }
    struct Cell { float a = 0, n = 0; };
    struct Entry { Cell c[NH]; };
    // размеры пространств ключей по семействам и смещения
    static constexpr int KS[NF] = {1, 41, 2, 82, 3, 123, 2, 3, 4, 5, 6, 7, 8, 41 * 41, 41 * 41 * 41};
    int off[NF + 1];
    std::vector<Entry> tab;
    float loss[NF][NH];
    float neff[NH];                           // затухающее число разрешённых наблюдений по горизонту
    int bestF[NH];                            // лучшее семейство по горизонту (пересчитывается раз в несколько разрешений)
    int sinceBest[NH];
    int curJ = -1;                            // выбранный адаптивно горизонт (-1: нет различающего признака)
    int needMask = (1 << NH) - 1;             // какие горизонты считаем
    int famMask = (1 << NF) - 1;              // какие семейства признаков используем
    // по трубам (pid с 1)
    std::vector<int> birth;
    std::vector<uint8_t> alive;
    std::vector<uint16_t> resolved;           // маска разрешённых горизонтов
    std::vector<int> kidx;                    // NF индексов записей на трубу, подряд
    size_t qhead[NH];                         // следующая труба (pid), которой исполняется возраст H_j
    int t = 0;                                // время = число обработанных команд
    int lp1 = 0, lp2 = 0, rp = 0;             // контекст: две прошлые длины, позиция в серии одинаковых длин
    int pm[9] = {0};                          // pid mod m для m = 2..8 (счётчики без деления)
    int strideT1 = 30000, strideT2 = 90000;   // границы времени для подвыборки обучающих труб

    Predictor() {
        off[0] = 0;
        for (int f = 0; f < NF; f++) off[f + 1] = off[f] + KS[f];
        tab.resize(off[NF]);
        for (int f = 0; f < NF; f++) for (int j = 0; j < NH; j++) loss[f][j] = 0.0f;
        for (int j = 0; j < NH; j++) { neff[j] = 0.0f; qhead[j] = 1; bestF[j] = 0; sinceBest[j] = 0; }
        birth.assign(1, 0); alive.assign(1, 0); resolved.assign(1, 0); kidx.assign(NF, 0);
    }

    // горизонт, на котором лучшее семейство признаков сильнее всего выигрывает у «общей» вероятности смерти
    int bestJ() const {
        int bj = -1;
        double bd = 0.02;   // минимальный выигрыш (нат на наблюдение), иначе предсказателя «нет»
        for (int j = 0; j < NH; j++) {
            if (!(needMask >> j & 1) || neff[j] < 15.0f) continue;
            float ml = loss[0][bestF[j]];
            double d = (double)(loss[0][j] - loss[bestF[j]][j]) / neff[j] * std::min(1.0, neff[j] / 40.0);
            (void)ml;
            if (d > bd) { bd = d; bj = j; }
        }
        return bj;
    }

    void recomputeBest(int j) {
        int bf = 0;
        float bl = 1e30f;
        for (int f = 0; f < NF; f++) if ((famMask >> f & 1) && loss[f][j] < bl) { bl = loss[f][j]; bf = f; }
        bestF[j] = bf;
        sinceBest[j] = 0;
    }

    // вызывается при укладке очередной трубы (до размещения); возвращает вероятности p[j] смерти за H_j
    void onPlace(int L, float* p) {
        if ((birth.size() & 15) == 0) curJ = bestJ();
        int pid = (int)birth.size();
        int rpn = (L == lp1) ? rp + 1 : 0;
        int r2 = rpn & 1, r3 = rpn % 3;
        for (int m = 2; m <= 8; m++) pm[m] = (pm[m] + 1 == m) ? 0 : pm[m] + 1;   // pm[m] = pid mod m
        int stride = t < strideT1 ? 1 : (t < strideT2 ? 2 : 4);   // на больших тестах таблицы обучаем по подвыборке труб
        birth.push_back(t); alive.push_back(1); resolved.push_back(stride > 1 && (pid % stride) != 0 ? (uint16_t)0xFFFF : (uint16_t)0);
        size_t base = kidx.size();
        kidx.resize(base + NF);
        int* k = &kidx[base];
        k[0] = off[0];
        k[1] = off[1] + L;
        k[2] = off[2] + r2;
        k[3] = off[3] + L * 2 + r2;
        k[4] = off[4] + r3;
        k[5] = off[5] + L * 3 + r3;
        for (int m = 2; m <= 8; m++) k[4 + m] = off[4 + m] + pm[m];
        k[13] = off[13] + lp1 * 41 + L;
        k[14] = off[14] + (lp2 * 41 + lp1) * 41 + L;
        for (int j = 0; j < NH; j++) {
            if (!(needMask >> j & 1)) { p[j] = 0.5f; continue; }
            const Cell& c = tab[k[bestF[j]]].c[j];
            p[j] = (c.a + 0.5f) / (c.n + 1.0f);
        }
        rp = rpn; lp2 = lp1; lp1 = L;
    }

    void resolve(int pid, int j, int y) {
        resolved[pid] |= (uint16_t)(1 << j);
        neff[j] = neff[j] * 0.995f + 1.0f;
        const int* k = &kidx[(size_t)pid * NF];
        for (int f = 0; f < NF; f++) {
            if (!(famMask >> f & 1)) continue;
            Cell& c = tab[k[f]].c[j];
            float p = (c.a + 0.5f) / (c.n + 1.0f);
            float l = -fastLog(y ? p : 1.0f - p);
            loss[f][j] = loss[f][j] * 0.995f + l;
            c.n += 1; c.a += (float)y;
        }
        if (++sinceBest[j] >= 8) recomputeBest(j);
    }

    void onDelete(int pid) {
        alive[pid] = 0;
        int age = t - birth[pid];
        for (int j = 0; j < NH; j++)
            if ((needMask >> j & 1) && !(resolved[pid] >> j & 1) && age <= H(j)) resolve(pid, j, 1);
    }

    // конец команды: продвинуть время, разрешить «дожила до H»
    void tick() {
        t++;
        size_t n = birth.size();
        for (int j = 0; j < NH; j++) {
            if (!(needMask >> j & 1)) continue;
            size_t& h = qhead[j];
            while (h < n && birth[h] + H(j) <= t) {
                int pid = (int)h++;
                if (!(resolved[pid] >> j & 1) && alive[pid]) resolve(pid, j, 0);
            }
        }
    }
};
constexpr int Predictor::KS[Predictor::NF];

// ---------------------------------------------------------------- движок (портфель)

struct EngineParams {   // параметры выбора лидера (стенд может менять)
    double hyst = 0.005;// переключаться, только если лидер лучше текущей на hyst (доля оценки)
    int absMargin = 0;  // ... и минимум на столько клеток
    int dwell = 500;    // минимум команд между переключениями
    double gamma = 0.5; // оценка = live + (1-gamma)*свободное + gamma*непригодное (дыры короче квантиля длин)
    int qwin = 256;     // окно (число последних укладок) для квантиля длин
    double qq = 0.25;   // квантиль
    int every = 8;      // оценивать лидера раз в столько укладок
    double hystArena = 0.03;  // гистерезис при переключении НА политику с аренами
    double hystAA = 0.03;     // ... если и текущая политика с аренами (переключение арена -> арена)
    int startCmds = 60;       // первые startCmds команд ведёт политика startName (хороша на крошечных «головоломках»), затем — политика 0
    double fast = 0.06;       // при выигрыше больше fast (доля оценки) переключаемся, не дожидаясь dwell
    int dwellFast = 40;       // ... но не чаще раза в dwellFast укладок
};
static EngineParams& engineParams() { static EngineParams ep; return ep; }
static bool& dbgSwitch() { static bool b = false; return b; }
static int g_famMask = -1, g_nhLimit = -1;   // стенд: ограничения предсказателя для замеров (в бою не используются)
static int& dbgScore() { static int v = 0; return v; }   // стенд: печатать оценки политик раз в dbgScore укладок
static bool& leanEnabled() { static bool b = true; return b; }   // стенд выключает (CPU-время процесса там общее для всех тестов)
#ifdef SIM
// Оракул времени смерти (только стенд): death[pid] = индекс команды удаления (или большое число)
static const std::vector<int>*& oracleDeath() { static thread_local const std::vector<int>* v = nullptr; return v; }
#endif

struct Engine {
    Alloc real;
    std::vector<Alloc> sh;   // теневые симуляции чистых политик (пусто в режиме одной политики)
    Predictor pred;
    int fixed;               // >= 0: одна чистая политика; -1: портфель
    int cur = 0;
    int sinceSwitch = 0;
    bool usePred = true;
    double lenHist[41] = {0};   // затухающая гистограмма длин (для таблицы W)
    double Wtab[WMAX + 1];      // W[s]: ожидаемый непригодный остаток в дыре размера s
    int ring[4096];          // длины последних укладок (окно квантиля)
    int ringCnt[41] = {0};
    int placed = 0;
    int placedAll = 0;
    int cmds = 0;
    int startIdx = 0;
    bool lean = false;       // аварийный режим: теней и предсказателя нет, реальная траншея ведётся одной политикой
    double leanPerCmd = LEAN_PER_CMD;   // порог CPU на команду: 5 с / 200000 = 25 мкс, берём запас

    explicit Engine(int p = -1) : fixed(p) {
        for (int i = 0; i <= WMAX; i++) Wtab[i] = i;
        for (int k = 0; k < NPOL; k++) if (!strcmp(POLS[k].name, "xw-R")) startIdx = k;
        if (p < 0 && engineParams().startCmds > 0) cur = startIdx;
        if (p < 0) sh.resize(NPOL);
        usePred = false;
        pred.needMask = 0;
        for (int k = 0; k < NPOL; k++) if ((p < 0 || k == p) && POLS[k].js >= 0) { usePred = true; pred.needMask |= POLS[k].js == 99 ? ((1 << Predictor::NH) - 1) : (1 << POLS[k].js); }
        if (g_famMask >= 0) pred.famMask = g_famMask | 1;
        if (g_nhLimit > 0) pred.needMask &= (1 << g_nhLimit) - 1;
    }
    int classOf(const Pol& p, const float* pr) const {
        if (p.js < 0) return 0;
        int j = p.js == 99 ? pred.curJ : p.js;
        return j >= 0 && pr[j] >= p.tau ? 1 : 0;
    }

    void updateW() {
        // W[s] = ожидаемый непригодный остаток, если дыру s заполнять трубами из наблюдаемого распределения (бесконечный запас)
        double tot = 0;
        for (int L = 1; L <= 40; L++) tot += lenHist[L];
        Wtab[0] = 0;
        for (int sI = 1; sI <= WMAX; sI++) {
            double z = 0, acc = 0;
            for (int L = 1; L <= 40 && L <= sI; L++) { z += lenHist[L]; acc += lenHist[L] * Wtab[sI - L]; }
            Wtab[sI] = (z <= 1e-9 * (tot + 1e-9) || tot <= 0) ? (double)sI : acc / z;
        }
    }

    double score(int k, int Lq) const {
        const EngineParams& ep = engineParams();
        const Alloc& a = sh[k];
        double freeK = (double)a.trench - (double)a.live;
        double unus = (double)(a.hs[0].unusable(Lq) + a.hs[1].unusable(Lq));
        return (double)a.live + (1.0 - ep.gamma) * freeK + ep.gamma * unus;
    }

    void goLean() {
        lean = true;
        int pick = 0;
        for (int k = 0; k < NPOL; k++) if (POLS[k].js < 0) { pick = k; if (k == cur) break; }
        if (POLS[cur].js < 0) pick = cur;
        fixed = pick;
        usePred = false;
        std::vector<Alloc>().swap(sh);
        pred = Predictor();
    }
    void place(int L, int& N, int& D) {
        if (leanEnabled() && !lean && (cmds & 4095) == 4095 && cmds >= 8192 && cpuSeconds() > leanPerCmd * cmds) goLean();
        float pr[Predictor::NH];
        for (int j = 0; j < Predictor::NH; j++) pr[j] = 0.5f;
        for (int l = 1; l <= 40; l++) lenHist[l] *= 0.999;
        lenHist[L] += 1.0;
        if ((placedAll++ & 127) == 0) updateW();
        real.W = Wtab;
        for (auto& a : sh) a.W = Wtab;
        if (usePred) pred.onPlace(L, pr);
#ifdef SIM
        if (oracleDeath()) {
            int pid = (int)real.pos.size();
            int life = (*oracleDeath())[pid] - cmds;
            for (int j = 0; j < Predictor::NH; j++) pr[j] = life <= Predictor::H(j) ? 1.0f : 0.0f;
        }
#endif
        cmds++;
        if (fixed >= 0) {
            real.place(POLS[fixed], L, classOf(POLS[fixed], pr), N, D);
        } else {
            const EngineParams& ep = engineParams();
            for (int k = 0; k < NPOL; k++) { int n, d; sh[k].place(POLS[k], L, classOf(POLS[k], pr), n, d); }
            // окно длин
            int w = std::min(ep.qwin, 4096);
            if (placed >= w) ringCnt[ring[placed % w]]--;
            ring[placed % w] = L; ringCnt[L]++; placed++;
            sinceSwitch++;
            if (placed == ep.startCmds && cur == startIdx && startIdx != 0) { cur = 0; sinceSwitch = 0; }
            if (placed % ep.every == 0 && sinceSwitch >= ep.dwellFast) {
                int n = std::min(placed, w), acc = 0, target = (int)(ep.qq * n), Lq;
                for (Lq = 1; Lq < 40; Lq++) { acc += ringCnt[Lq]; if (acc > target) break; }
                int best = cur;
                double sBest = score(cur, Lq), sCur = sBest;
                for (int k = 0; k < NPOL; k++) { double sk = score(k, Lq); if (sk < sBest) { sBest = sk; best = k; } }
                if (dbgScore() && placed % dbgScore() == 0) {
                    fprintf(stderr, "t=%d cur=%s Lq=%d scores:", cmds, POLS[cur].name, Lq);
                    for (int k = 0; k < NPOL; k++) fprintf(stderr, " %.0f", score(k, Lq));
                    fputc(10, stderr);
                }
                double lead = sCur > 0 ? (sCur - sBest) / sCur : 0;
                double need = POLS[best].js >= 0 ? (POLS[cur].js >= 0 ? ep.hystAA : ep.hystArena) : ep.hyst;
                if (best != cur && lead > need && (sinceSwitch >= ep.dwell || lead > ep.fast) && sCur - sBest > ep.absMargin) {
                    if (dbgSwitch()) {
                        fprintf(stderr, "t=%d switch %s -> %s  trench: real=%d", cmds, POLS[cur].name, POLS[best].name, real.trench);
                        for (int k = 0; k < NPOL; k++) fprintf(stderr, " %s=%d", POLS[k].name, sh[k].trench);
                        fputc('\n', stderr);
                    }
                    cur = best;
                    sinceSwitch = 0;
                }
            }
            real.place(POLS[cur], L, classOf(POLS[cur], pr), N, D);
        }
        if (usePred) pred.tick();
    }
    int remove(int pid) {
        cmds++;
        if (usePred) { pred.onDelete(pid); pred.tick(); }
        for (auto& a : sh) a.remove(pid);
        return real.remove(pid);
    }
};

#ifndef SIM
// ---------------------------------------------------------------- ввод-вывод (боевой режим)
#if defined(_MSC_VER)
#define GETC _getchar_nolock
#elif defined(_WIN32)
#define GETC getchar
#else
#define GETC getchar_unlocked
#endif
static bool readInt(int& x) {
    int c = GETC();
    while (c != EOF && c != '-' && (c < '0' || c > '9')) c = GETC();
    if (c == EOF) return false;
    bool neg = false;
    if (c == '-') { neg = true; c = GETC(); }
    int v = 0;
    while (c >= '0' && c <= '9') { v = v * 10 + (c - '0'); c = GETC(); }
    x = neg ? -v : v;
    return true;
}
static void writeInts(int a, int b, bool two) {
    char buf[32];
    int n = 0;
    auto put = [&](int v) {
        char t[12]; int k = 0;
        if (v == 0) t[k++] = '0';
        while (v > 0) { t[k++] = (char)('0' + v % 10); v /= 10; }
        while (k) buf[n++] = t[--k];
    };
    put(a);
    if (two) { buf[n++] = ' '; put(b); }
    buf[n++] = '\n';
    fwrite(buf, 1, n, stdout);
    fflush(stdout);
}
#ifndef FIXED_POL
#define FIXED_POL -1
#endif
int main() {
    Engine eng(FIXED_POL);
    int c;
    while (readInt(c)) {
        if (c == 0) return 0;
        if (c > 0) {
            int N, D;
            eng.place(c, N, D);
            writeInts(N, D, true);
        } else {
            writeInts(eng.remove(-c), 0, false);
        }
    }
}
#endif
