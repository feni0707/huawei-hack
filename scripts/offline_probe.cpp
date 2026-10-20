// Оффлайн-зонд: проверка ДОСТИЖИМОСТИ тривиальной верхней оценки UB = live_end / (20*ceil(peak/20)).
//
// Оффлайн-постановка: будущее известно, копать можно заранее в любом объёме, траншея не укорачивается =>
// задача сводится ровно к упаковке интервалов: пересекающиеся по времени жизни трубы не пересекаются по клеткам,
// минимизируется max(pos + len). Это classic dynamic storage allocation (NP-трудная в общем виде),
// поэтому проверяем двумя жадными оффлайн-алгоритмами, знающими время смерти каждой трубы:
//   first = left edge (первое подходящее место слева), best = самая тесная подходящая дыра.
// Если оффлайн-жадина упирается в UB — значит UB достижима оффлайном.
//
//   build: cl /nologo /O2 /EHsc /std:c++17 /Fe:build\offline_probe.exe scripts\offline_probe.cpp
//   run:   build\offline_probe.exe materials\tests\*.in
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <vector>

typedef long long ll;

static std::vector<ll> readTest(const char* fn, double& w) {
    FILE* f = fopen(fn, "r");
    if (!f) { perror(fn); exit(1); }
    int n;
    if (fscanf(f, "%lf %d", &w, &n) != 2) exit(1);
    std::vector<ll> c(n);
    for (int i = 0; i < n; i++) if (fscanf(f, "%lld", &c[i]) != 1) exit(1);
    fclose(f);
    return c;
}

static const int MAXLEN = 40;

// Одна реализация: раскладка по времени с известными смертями.
struct Packer {
    int mode;                                    // 0 = first fit (left edge), 1 = best fit
    std::set<std::pair<ll, ll>> gapsByPos;       // (start, len)
    std::set<std::pair<ll, ll>> gapsByLen;       // (len, start)
    std::map<ll, ll> gapByEnd;                   // end -> start  (для слияния при освобождении)
    std::multimap<int, ll> deaths;               // death time -> start размещённого блока
    std::vector<std::set<ll>> avail;             // avail[L] = начала дыр длины >= L (L = 1..40)
    ll used = 0;                                 // верхняя граница занятой области

    void init(int m) {
        mode = m;
        gapsByPos.clear(); gapsByLen.clear(); gapByEnd.clear(); deaths.clear();
        avail.assign(MAXLEN + 1, {});
        used = 0;
    }

    void addGap(ll s, ll len) {
        if (len <= 0) return;
        gapsByPos.insert({ s, len });
        gapsByLen.insert({ len, s });
        gapByEnd[s + len] = s;
        for (int L = 1; L <= MAXLEN && L <= len; L++) avail[L].insert(s);
    }
    void delGap(ll s, ll len) {
        gapsByPos.erase({ s, len });
        gapsByLen.erase({ len, s });
        gapByEnd.erase(s + len);
        for (int L = 1; L <= MAXLEN && L <= len; L++) avail[L].erase(s);
    }

    ll findGap(int L) {
        if (mode == 0) {
            return avail[L].empty() ? -1 : *avail[L].begin();
        }
        auto it = gapsByLen.lower_bound({ L, -1 });
        return it == gapsByLen.end() ? -1 : it->second;
    }

    ll lenOfGapStart(ll s) {
        auto it = gapsByPos.lower_bound({ s, -1 });
        return it == gapsByPos.end() || it->first != s ? -1 : it->second;
    }

    void place(int L, int deathAt, int now, std::vector<ll>& startOf);

    // освободить блок [p, p+len): слить с соседними дырами
    void release(ll p, int len) {
        ll s = p, e = p + (ll)len;
        auto itEnd = gapByEnd.find(s);
        if (itEnd != gapByEnd.end()) {
            ll gs = itEnd->second;
            delGap(gs, s - gs);
            s = gs;
        }
        auto itPos = gapsByPos.lower_bound({ e, -1 });
        if (itPos != gapsByPos.end() && itPos->first == e) {
            ll gl = itPos->second;
            delGap(e, gl);
            e += gl;
        }
        addGap(s, e - s);
    }
};

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s tests/*.in\n", argv[0]); return 1; }
    double tot[2] = { 0, 0 }, totUB = 0, totW = 0;
    printf("%-30s %4s %9s %9s %9s %9s %9s  %s\n", "test", "w", "live_end", "peak", "Tmin(UB)",
           "T_first", "T_best", "reach");
    for (int a = 1; a < argc; a++) {
        double w;
        std::vector<ll> cmd = readTest(argv[a], w);
        int n = (int)cmd.size();
        std::vector<int> plen(1, 0), death(1, n);
        int pid = 0;
        for (int i = 0; i < n; i++) {
            if (cmd[i] > 0) { pid++; plen.push_back((int)cmd[i]); death.push_back(n); }
            else death[(int)-cmd[i]] = i;
        }
        ll live = 0, peak = 0;
        for (int i = 0; i < n; i++) {
            if (cmd[i] > 0) { live += cmd[i]; peak = std::max(peak, live); }
            else live -= plen[(int)-cmd[i]];
        }
        ll Tmin = 20 * ((peak + 19) / 20);

        ll T[2] = { 0, 0 };
        for (int mode = 0; mode < 2; mode++) {
            Packer P;
            P.init(mode);
            std::vector<ll> pos(pid + 1, -1);
            int curPid = 0;
            // активные трубы: death time -> (pos, len)
            std::multimap<int, std::pair<ll, int>> dying;
            for (int i = 0; i < n; i++) {
                while (!dying.empty() && dying.begin()->first <= i) {
                    auto pr = dying.begin()->second;
                    dying.erase(dying.begin());
                    P.release(pr.first, pr.second);
                }
                if (cmd[i] <= 0) continue;
                int L = (int)cmd[i];
                int id = ++curPid;
                ll p = P.findGap(L);
                if (p >= 0) {
                    ll len = P.lenOfGapStart(p);
                    P.delGap(p, len);
                    if (len > L) P.addGap(p + L, len - L);
                } else {
                    p = P.used;
                    P.used += L;
                }
                P.used = std::max(P.used, p + L);
                pos[id] = p;
                if (death[id] < n) dying.insert({ death[id], { p, L } });
            }
            T[mode] = P.used;
        }
        ll Tf = 20 * ((T[0] + 19) / 20), Tb = 20 * ((T[1] + 19) / 20);
        double ub = (double)live / Tmin, kf = (double)live / Tf, kb = (double)live / Tb;
        totUB += ub * w; tot[0] += kf * w; tot[1] += kb * w; totW += w;
        printf("%-30s %4.1f %9lld %9lld %9lld %9lld %9lld  %s\n", argv[a], w, live, peak, Tmin, Tf, Tb,
               (Tf <= Tmin ? "FIRST=UB" : (Tb <= Tmin ? "BEST=UB" : "-")));
    }
    printf("\nweighted of %.1f: UB=%.4f  first=%.4f  best=%.4f\n", totW, totUB, tot[0], tot[1]);
    return 0;
}
