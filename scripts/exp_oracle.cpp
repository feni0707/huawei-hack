// Эксперимент: насколько помогает знание времени смерти труб (оракул) при размещении по «похожести соседей».
// Это потолок для любого предсказателя времени жизни. Полный перебор дыр (медленно, только для оценки).
//   exp_oracle <wm> <wr> <wt> <test.in>...
#define SIM
#include "../solutions/portfolio.cpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>

struct XAlloc : Alloc {
    std::map<int, int> pipeAt;  // начало -> pid (живые трубы)

    void putExplicit(int L, int D, int N) {
        trench += 20 * N;
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
        int pid = (int)pos.size();
        pos.push_back(D);
        len.push_back(L);
        live += L;
        pipeAt[D] = pid;
    }
    int removeX(int pid) {
        pipeAt.erase(pos[pid]);
        return remove(pid);
    }
    // pid трубы, заканчивающейся ровно в x (левый сосед), или 0
    int leftNeighbor(int x) const {
        if (x <= 0) return 0;
        auto it = pipeAt.lower_bound(x);
        if (it == pipeAt.begin()) return 0;
        --it;
        return it->second;  // проверка вплотную — снаружи
    }
    int rightNeighbor(int x) const {
        auto it = pipeAt.find(x);
        return it == pipeAt.end() ? 0 : it->second;
    }
};

int main(int argc, char** argv) {
    double wm = atof(argv[1]), wr = atof(argv[2]), wt = atof(argv[3]);
    double tot = 0;
    for (int ai = 4; ai < argc; ai++) {
        std::ifstream f(argv[ai]);
        double w; int n;
        f >> w >> n;
        std::vector<int> cmds(n);
        for (auto& c : cmds) f >> c;
        int np = 0;
        for (int c : cmds) if (c > 0) np++;
        const int INF = 1 << 28;
        std::vector<int> death(np + 2, INF);
        for (int i = 0; i < n; i++) if (cmds[i] < 0) death[-cmds[i]] = i;
        XAlloc A;
        auto remaining = [&](int pid, int now) {  // log2(1+остаток жизни)
            if (pid == 0) return 0.0;
            int d = death[pid];
            if (d >= INF) return 24.0;
            return std::log2(1.0 + std::max(0, d - now));
        };
        for (int i = 0; i < n; i++) {
            int c = cmds[i];
            if (c < 0) { A.removeX(-c); continue; }
            int pid = (int)A.pos.size();
            double up = remaining(pid, i);
            int tailLen = A.trench - A.tailStart;
            // кандидаты: (стоимость, D, N)
            double bestCost = 1e18; int bestD = -1, bestN = 0;
            auto mism = [&](int nb) { return nb ? std::fabs(up - remaining(nb, i)) : 0.0; };
            auto consider = [&](int D, int N, int h /*размер дыры или -1 для хвоста*/, bool leftSide) {
                double cost = 0;
                int adj;
                if (h >= 0) {
                    int a = leftSide ? D : D - (h - c);  // начало дыры
                    int r = h - c;
                    if (r == 0) cost -= 1.0 * wr * 100;
                    else cost += wr * (r < 8 ? 1.0 : 0.0) * (8 - std::min(r, 8));
                    adj = leftSide ? A.leftNeighbor(a) : A.rightNeighbor(a + h);
                } else {
                    cost += wt * (20.0 * N);
                    adj = A.leftNeighbor(D);
                }
                cost += wm * mism(adj);
                if (cost < bestCost) { bestCost = cost; bestD = D; bestN = N; }
            };
            for (auto& kv : A.byAddr) {
                int a = kv.first, h = kv.second;
                if (h < c) continue;
                consider(a, 0, h, true);
                if (h > c) consider(a + h - c, 0, h, false);
            }
            if (tailLen >= c) consider(A.tailStart, 0, -1, true);
            else consider(A.tailStart, (c - tailLen + 19) / 20, -1, true);
            A.putExplicit(c, bestD, bestN);
        }
        double kpi = (double)A.live / A.trench;
        printf("%s kpi=%.4f trench=%d live=%lld\n", argv[ai], kpi, A.trench, A.live);
        tot += kpi * w;
    }
    printf("sum w*kpi = %.4f\n", tot);
}
