// Эксперимент: арены (чанки) по классам времени жизни. Класс — по истинной смерти (оракул) или по предсказателю.
//   exp_arena <H> <C> <mode> <test.in>...   mode: 0 = оракул (смерть ближе H команд → transient), 1 = bf для сравнения
#define SIM
#include "../solutions/portfolio.cpp"

#include <fstream>
#include <iostream>

struct Arena {
    int x = 0, y = 0;   // свободный диапазон [x,y), зарезервированный классом (пуст, если x==y)
};

struct CAlloc : Alloc {
    Arena ar[2];
    int C = 100;

    // положить в явную позицию D (аналог placeExplicit)
    void putAt(int L, int D) {
        // D лежит либо в дыре, либо в хвосте, либо в зарезервированной арене (арена вне byAddr/хвоста)
        if (D >= tailStart) {
            if (D > tailStart) addHole(tailStart, D - tailStart);
            tailStart = D + L;
        } else {
            auto it = byAddr.upper_bound(D);
            --it;
            int a = it->first, h = it->second;
            delHole(a, h);
            if (D > a) addHole(a, D - a);
            if (D + L < a + h) addHole(D + L, a + h - (D + L));
        }
        pos.push_back(D);
        len.push_back(L);
        live += L;
    }
    // зарезервировать диапазон [a, a+size) из дыры (a,h) целиком или частично — вырезаем из пула
    void reserve(int a, int size) {
        auto it = byAddr.upper_bound(a);
        --it;
        int ha = it->first, h = it->second;
        delHole(ha, h);
        if (a > ha) addHole(ha, a - ha);
        if (a + size < ha + h) addHole(a + size, ha + h - (a + size));
    }
    void reserveTail(int a, int size) {  // [a, a+size) в хвосте (a >= tailStart)
        if (a > tailStart) addHole(tailStart, a - tailStart);
        tailStart = a + size;
    }
    void release(int cls) {
        Arena& A = ar[cls];
        if (A.x < A.y) {
            // вернуть остаток арены в пул: аналог remove-слияния; делаем через временную «трубу»
            int s = A.x, e = A.y;
            auto it = byAddr.lower_bound(s);
            if (it != byAddr.begin()) {
                auto jt = std::prev(it);
                if (jt->first + jt->second == s) { s = jt->first; int hh = jt->second; delHole(jt->first, hh); }
            }
            if (e == tailStart) {
                tailStart = s;
            } else {
                auto r = byAddr.find(e);
                if (r != byAddr.end()) { int hh = r->second; e += hh; delHole(r->first, hh); }
                addHole(s, e - s);
            }
        }
        A.x = A.y = 0;
    }

    // N докопки возвращается через trench-изменение
    void place(int L, int cls) {
        Arena& A = ar[cls];
        if (A.y - A.x >= L) {
            int D = A.x;
            A.x += L;
            pos.push_back(D); len.push_back(L); live += L;
            return;
        }
        release(cls);
        // 1) скрап: маленькие дыры (< C) — best-fit, но только если подходят
        // 2) новая арена: самая большая дыра >= max(C, L); иначе хвост
        int a = -1, h = 0;
        if (!big.empty() || mask >> 40 & 1) {
            a = big.rbegin()->second; h = big.rbegin()->first;
        }
        // для persistent (cls=1): сначала точный best-fit в скрапе
        if (cls == 1 || true) {
            uint64_t mm = mask & (~0ull << L) & ((1ull << 40) - 1);
            // скрап подходит, если размер < C
            int s = mm ? ctz64(mm) : -1;
            if (s >= 0 && s < C && s < 40) {
                int aa = *sz[s].begin();
                putAt(L, aa);
                return;
            }
        }
        int need = std::max(C, L);
        if (a >= 0 && h >= need) {
            int size = std::min(h, std::max(need, C));
            reserve(a, size);
            A.x = a; A.y = a + size;
        } else {
            int tailLen = trench - tailStart;
            if (tailLen < need) trench += 20 * ((need - tailLen + 19) / 20);
            reserveTail(tailStart, need);
            A.x = tailStart - need; A.y = A.x + need;
        }
        int D = A.x;
        A.x += L;
        pos.push_back(D); len.push_back(L); live += L;
    }
    int removeP(int pid) { return remove(pid); }
};

int main(int argc, char** argv) {
    int H = atoi(argv[1]), C = atoi(argv[2]), mode = atoi(argv[3]);
    double tot = 0;
    for (int ai = 4; ai < argc; ai++) {
        std::ifstream f(argv[ai]);
        double w; int n;
        f >> w >> n;
        std::vector<int> cmds(n);
        for (auto& c : cmds) f >> c;
        int np = 0;
        for (int c : cmds) if (c > 0) np++;
        std::vector<int> death(np + 2, 1 << 28);
        for (int i = 0; i < n; i++) if (cmds[i] < 0) death[-cmds[i]] = i;
        CAlloc A; A.C = C;
        Alloc B;
        for (int i = 0; i < n; i++) {
            int c = cmds[i];
            if (mode == 1) {
                if (c > 0) { int N, D; B.place(POLS[0], c, N, D); } else B.remove(-c);
                continue;
            }
            if (c < 0) { A.removeP(-c); continue; }
            int pid = (int)A.pos.size();
            int cls = (death[pid] - i <= H) ? 0 : 1;
            A.place(c, cls);
        }
        double kpi = mode == 1 ? (double)B.live / B.trench : (double)A.live / A.trench;
        // арены, не возвращённые, — потери, но они в траншее (учтено)
        printf("%s kpi=%.4f\n", argv[ai], kpi);
        tot += kpi * w;
    }
    printf("sum w*kpi = %.4f\n", tot);
}
