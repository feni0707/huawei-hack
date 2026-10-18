// Эксперимент 2: дыры с тегами. T-класс (скоротечные) размещается best-fit по T-дырам, иначе открывает «чанк» в хвосте;
// P-класс (долгоживущие) — best-fit по U-дырам и малым T-дырам, иначе хвост. Класс — по оракулу смерти.
//   exp_tag <H> <C> <Tsmall> <mode> <test.in>...
#define SIM
#include "../solutions/portfolio.cpp"

#include <fstream>
#include <iostream>

struct HoleSet {
    std::set<int> sz[40];
    std::set<std::pair<int, int>> big;
    uint64_t mask = 0;
    void add(int a, int h) {
        if (h < 40) { sz[h].insert(a); mask |= 1ull << h; }
        else { big.insert({h, a}); mask |= 1ull << 40; }
    }
    void del(int a, int h) {
        if (h < 40) { sz[h].erase(a); if (sz[h].empty()) mask &= ~(1ull << h); }
        else { big.erase({h, a}); if (big.empty()) mask &= ~(1ull << 40); }
    }
    // best fit: наименьший размер >= L, меньший адрес; maxSize ограничивает размер дыры (для P в T-дырах)
    bool bestFit(int L, int maxSize, int& a, int& h) const {
        uint64_t mm = mask & (~0ull << L);
        if (!mm) return false;
        int s = ctz64(mm);
        if (s < 40) { if (s > maxSize) return false; a = *sz[s].begin(); h = s; return true; }
        if (big.begin()->first > maxSize) return false;
        a = big.begin()->second; h = big.begin()->first; return true;
    }
};

struct TAlloc {
    int trench = 0, tailStart = 0;
    std::map<int, std::pair<int, int>> byAddr;  // начало -> (длина, тег)
    HoleSet hs[2];                              // 0 = U, 1 = T
    std::vector<int> pos, len, cls;
    long long live = 0;
    int C = 100, Tsmall = 40, Uuse = 0;
    TAlloc() { pos.assign(1, -1); len.assign(1, 0); cls.assign(1, 0); }

    void addHole(int a, int h, int tag) { byAddr[a] = {h, tag}; hs[tag].add(a, h); }
    void delHole(int a) { auto it = byAddr.find(a); hs[it->second.second].del(a, it->second.first); byAddr.erase(it); }

    void carve(int a, int h, int tag, int D, int L) {  // из дыры (a,h,tag) вырезать [D, D+L)
        delHole(a);
        if (D > a) addHole(a, D - a, tag);
        if (D + L < a + h) addHole(D + L, a + h - (D + L), tag);
    }
    void commit(int D, int L, int c) { pos.push_back(D); len.push_back(L); cls.push_back(c); live += L; }
    int tailPlace(int L, int c, int& N) {
        int tailLen = trench - tailStart;
        N = 0;
        if (tailLen < L) { N = (L - tailLen + 19) / 20; trench += 20 * N; }
        int D = tailStart;
        tailStart += L;
        commit(D, L, c);
        return D;
    }

    void place(int L, int c, int& N, int& D) {
        N = 0;
        int a, h;
        if (c == 1) {
            if (hs[1].bestFit(L, 1 << 30, a, h)) { D = a; carve(a, h, 1, D, L); commit(D, L, c); return; }
            if (Uuse && hs[0].bestFit(L, L + Uuse, a, h)) { D = a; carve(a, h, 0, D, L); commit(D, L, c); return; }
            // новый чанк в хвосте
            int need = std::max(C, L);
            int tailLen = trench - tailStart;
            if (tailLen < need) { N = (need - tailLen + 19) / 20; trench += 20 * N; }
            int cs = tailStart;
            tailStart += need;
            if (need > L) addHole(cs + L, need - L, 1);
            D = cs;
            commit(D, L, c);
            return;
        }
        // P: best-fit по U-дырам и малым T-дырам
        int a0, h0, a1, h1;
        bool f0 = hs[0].bestFit(L, 1 << 30, a0, h0), f1 = hs[1].bestFit(L, Tsmall, a1, h1);
        if (f0 || f1) {
            bool useT = f1 && (!f0 || h1 < h0);
            if (useT) { D = a1; carve(a1, h1, 1, D, L); } else { D = a0; carve(a0, h0, 0, D, L); }
            commit(D, L, c);
            return;
        }
        D = tailPlace(L, c, N);
    }

    void remove(int pid) {
        int p = pos[pid], l = len[pid], tag = cls[pid];
        pos[pid] = -1;
        live -= l;
        int s = p, e = p + l;
        auto it = byAddr.lower_bound(p);
        if (it != byAddr.begin()) {
            auto jt = std::prev(it);
            if (jt->first + jt->second.first == p) { s = jt->first; tag = std::max(tag, jt->second.second); delHole(jt->first); }
        }
        if (e == tailStart) { tailStart = s; return; }
        auto r = byAddr.find(e);
        if (r != byAddr.end()) { tag = std::max(tag, r->second.second); e += r->second.first; delHole(r->first); }
        addHole(s, e - s, tag);
    }
};

int main(int argc, char** argv) {
    int H = atoi(argv[1]), C = atoi(argv[2]), Ts = atoi(argv[3]), Uuse = atoi(argv[4]);
    double tot = 0;
    for (int ai = 5; ai < argc; ai++) {
        std::ifstream f(argv[ai]);
        double w; int n;
        f >> w >> n;
        std::vector<int> cmds(n);
        for (auto& c : cmds) f >> c;
        int np = 0;
        for (int c : cmds) if (c > 0) np++;
        std::vector<int> death(np + 2, 1 << 28);
        for (int i = 0; i < n; i++) if (cmds[i] < 0) death[-cmds[i]] = i;
        TAlloc A; A.C = C; A.Tsmall = Ts; A.Uuse = Uuse;
        for (int i = 0; i < n; i++) {
            int c = cmds[i];
            if (c < 0) { A.remove(-c); continue; }
            int pid = (int)A.pos.size();
            int cls = (death[pid] - i <= H) ? 1 : 0;
            int N, D;
            A.place(c, cls, N, D);
        }
        double kpi = (double)A.live / A.trench;
        printf("%s=%.4f ", argv[ai], kpi);
        tot += kpi * w;
    }
    printf("\nsum w*kpi = %.4f\n", tot);
}
