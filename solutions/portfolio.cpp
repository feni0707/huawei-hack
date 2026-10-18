// Вариант B: портфель политик размещения + выбор политики по истории (теневые симуляции).
// Сборка решения:  c++ -std=c++17 -O2 -o build/portfolio solutions/portfolio.cpp
// Стенд (внутрипроцессный прогон тестов): scripts/sim.cpp, сборка scripts/build_sim.sh
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <vector>
#ifdef _MSC_VER
#include <intrin.h>
#endif

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

// ---------------------------------------------------------------- политики

enum Sel { BEST, FIRST, LAST, WORST, EXACTWORST };

struct Pol {
    const char* name;
    Sel sel;
    bool tieHigh;    // при равенстве размеров брать дыру с большим адресом
    bool right;      // класть к правому краю дыры
    bool tailFirst;  // хвост — обычный кандидат: если влезает, класть туда до поиска дыр
};

static const Pol POLS[] = {
    {"bf-L", BEST, false, false, false},    {"bf-R", BEST, false, true, false},
    {"bf-hiL", BEST, true, false, false},   {"bf-hiR", BEST, true, true, false},
    {"ff-L", FIRST, false, false, false},   {"ff-R", FIRST, false, true, false},
    {"lf-L", LAST, false, false, false},    {"lf-R", LAST, false, true, false},
    {"wf-L", WORST, false, false, false},   {"xw-L", EXACTWORST, false, false, false},
    {"bf-Lt", BEST, false, false, true},    {"ff-Lt", FIRST, false, false, true},
};
static const int NPOL = sizeof(POLS) / sizeof(POLS[0]);

// ---------------------------------------------------------------- состояние траншеи

struct Alloc {
    int trench = 0;      // длина траншеи
    int tailStart = 0;   // начало свободного хвоста (хвост = [tailStart, trench))
    std::map<int, int> byAddr;            // свободные дыры внутри траншеи (не касаются конца): начало -> длина
    std::set<int> sz[40];                 // дыры длины 1..39: адреса
    std::set<std::pair<int, int>> big;    // дыры длины >= 40: (длина, адрес)
    std::set<int> bigAddr;                // они же по адресу
    uint64_t mask = 0;                    // бит s (1..39) — есть дыра длины s; бит 40 — есть большая
    std::vector<int> pos, len;            // по PID (с 1)
    long long live = 0;

    Alloc() { pos.assign(1, -1); len.assign(1, 0); }

    void addHole(int a, int h) {
        byAddr[a] = h;
        if (h < 40) { sz[h].insert(a); mask |= 1ull << h; }
        else { big.insert({h, a}); bigAddr.insert(a); mask |= 1ull << 40; }
    }
    void delHole(int a, int h) {
        byAddr.erase(a);
        if (h < 40) { sz[h].erase(a); if (sz[h].empty()) mask &= ~(1ull << h); }
        else { big.erase({h, a}); bigAddr.erase(a); if (big.empty()) mask &= ~(1ull << 40); }
    }

    // кандидат-дыра по правилу политики; false — ни одна не подходит
    bool pick(const Pol& p, int L, int& a, int& h) const {
        uint64_t mm = mask & (~0ull << L);
        if (!mm) return false;
        switch (p.sel) {
        case BEST: {
            int s = ctz64(mm);
            if (s < 40) { a = p.tieHigh ? *sz[s].rbegin() : *sz[s].begin(); h = s; }
            else {
                auto it = big.begin();
                if (p.tieHigh) { it = big.upper_bound({it->first, 0x7fffffff}); --it; }
                a = it->second; h = it->first;
            }
            return true;
        }
        case FIRST:
        case LAST: {
            bool first = p.sel == FIRST;
            int best = first ? 0x7fffffff : -1, bh = 0;
            for (uint64_t m = mm & ((1ull << 40) - 1); m; m &= m - 1) {
                int s = ctz64(m);
                int x = first ? *sz[s].begin() : *sz[s].rbegin();
                if (first ? x < best : x > best) { best = x; bh = s; }
            }
            if (mm >> 40 & 1) {
                int x = first ? *bigAddr.begin() : *bigAddr.rbegin();
                if (first ? x < best : x > best) { best = x; bh = byAddr.find(x)->second; }
            }
            a = best; h = bh;
            return true;
        }
        case WORST:
            if (mm >> 40 & 1) { a = big.rbegin()->second; h = big.rbegin()->first; }
            else { int s = clz64idx(mm); a = *sz[s].begin(); h = s; }
            return true;
        case EXACTWORST:
            if (L < 40 && (mask >> L & 1)) { a = *sz[L].begin(); h = L; return true; }
            if (mm >> 40 & 1) { a = big.rbegin()->second; h = big.rbegin()->first; }
            else { int s = clz64idx(mm); a = *sz[s].begin(); h = s; }
            return true;
        }
        return false;
    }

    // уложить трубу длины L; возвращает N (докопка в 20-метровых сегментах) и D (позиция)
    void place(const Pol& p, int L, int& N, int& D) {
        N = 0;
        int a, h;
        int tailLen = trench - tailStart;
        if (p.tailFirst && tailLen >= L) {
            D = tailStart;
            tailStart += L;
        } else if (pick(p, L, a, h)) {
            delHole(a, h);
            if (p.right) { D = a + h - L; if (h > L) addHole(a, h - L); }
            else { D = a; if (h > L) addHole(a + L, h - L); }
        } else {
            if (tailLen < L) { N = (L - tailLen + 19) / 20; trench += 20 * N; }
            D = tailStart;
            tailStart += L;
        }
        int pid = (int)pos.size();
        pos.push_back(D);
        len.push_back(L);
        live += L;
    }

    int remove(int pid) {
        int p = pos[pid], l = len[pid];
        pos[pid] = -1;
        live -= l;
        int s = p, e = p + l;
        auto it = byAddr.lower_bound(p);
        if (it != byAddr.begin()) {
            auto jt = std::prev(it);
            if (jt->first + jt->second == p) { s = jt->first; int hh = jt->second; delHole(jt->first, hh); }
        }
        if (e == tailStart) {
            tailStart = s;
        } else {
            auto r = byAddr.find(e);
            if (r != byAddr.end()) { int hh = r->second; e += hh; delHole(r->first, hh); }
            addHole(s, e - s);
        }
        return p;
    }
};

// ---------------------------------------------------------------- движок (портфель)

struct Engine {
    Alloc real;
    std::vector<Alloc> sh;   // теневые симуляции чистых политик (пусто в режиме одной политики)
    int fixed;               // >= 0: одна чистая политика; -1: портфель
    int cur = 0;
    double hyst = 0.0;       // гистерезис: переключаться, только если лидер короче текущей на hyst*trench

    explicit Engine(int p = -1) : fixed(p) {
        if (p < 0) sh.resize(NPOL);
    }
    void place(int L, int& N, int& D) {
        if (fixed >= 0) { real.place(POLS[fixed], L, N, D); return; }
        for (int k = 0; k < NPOL; k++) { int n, d; sh[k].place(POLS[k], L, n, d); }
        int best = cur;
        for (int k = 0; k < NPOL; k++) if (sh[k].trench < sh[best].trench) best = k;
        if (best != cur && sh[best].trench < sh[cur].trench - hyst * sh[cur].trench) cur = best;
        real.place(POLS[cur], L, N, D);
    }
    int remove(int pid) {
        for (auto& a : sh) a.remove(pid);
        return real.remove(pid);
    }
};

#ifndef SIM
// ---------------------------------------------------------------- ввод-вывод (боевой режим)
int main() {
    Engine eng(-1);
    int c;
    while (scanf("%d", &c) == 1) {
        if (c == 0) return 0;
        if (c > 0) {
            int N, D;
            eng.place(c, N, D);
            printf("%d %d\n", N, D);
        } else {
            printf("%d\n", eng.remove(-c));
        }
        fflush(stdout);
    }
}
#endif
