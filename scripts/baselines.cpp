// Простые эталонные стратегии для стенда scripts/simc.cpp (Alloc-совместимый класс), режим — env K="0=<mode>":
//   0 — best fit (самый короткий промежуток ≥ L, среди равных левый; хвост — только если внутри не влезает)
//   1 — first fit (самый левый промежуток ≥ L)
//   2 — оракул-BF: среди промежутков ≥ L — best fit, но точное попадание и остаток кладутся к соседу, который
//       умрёт раньше (знание будущего через ORD — время удаления по PID), «вечные» трубы — в точные дыры.
//   c++ -std=c++17 -O2 -DSFSRC='"'$PWD/scripts/baselines.cpp'"' -o build/simc_base scripts/simc.cpp
#include <map>
#include <set>
#include <vector>
typedef long long ll;
extern double K[64];
extern std::vector<int> ORD; extern long long CUR;
struct ZStub { int get(ll) const { return 0; } };
struct Alloc {
    ll T = 0;
    ZStub zl;
    std::map<ll, ll> gs;              // start -> len (все свободные промежутки, включая хвост)
    std::set<std::pair<ll, ll>> byLen; // (len, start)
    std::vector<ll> pos; std::vector<int> len, pidAt;
    std::vector<int> owner;
    Alloc() { pos.push_back(-1); len.push_back(0); }
    void add(ll s, ll l) { gs[s] = l; byLen.insert({l, s}); }
    void del(ll s) { auto it = gs.find(s); byLen.erase({it->second, s}); gs.erase(it); }
    void freeR(ll s, ll l) {
        ll e = s + l;
        auto it = gs.find(e);
        if (it != gs.end()) { ll r = it->second; del(e); e += r; }
        it = gs.lower_bound(s);
        if (it != gs.begin()) { --it; if (it->first + it->second == s) { s = it->first; del(s); } }
        add(s, e - s);
    }
    void take(ll g, ll d, ll l) {
        ll gl = gs[g]; del(g);
        if (d > g) add(g, d - g);
        if (d + l < g + gl) add(d + l, g + gl - d - l);
    }
    ll tailStart() { if (gs.empty()) return T; auto it = std::prev(gs.end()); return it->first + it->second == T ? it->first : T; }
    int nb(ll x) { return (x >= 0 && x < T) ? owner[x] : 0; }
    ll death(int p) { return p > 0 && p < (int)ORD.size() ? ORD[p] : -1; }  // -1 — край
    void place(int L, ll& N, ll& D) {
        int mode = (int)K[0];
        ll ts = tailStart();
        ll bg = -1, bl = 0;
        if (mode == 1) {
            for (auto& [s, l] : gs) if (s != ts && l >= L) { bg = s; bl = l; break; }
        } else if (mode == 0) {
            for (auto it = byLen.lower_bound({L, -1}); it != byLen.end(); ++it) if (it->second != ts) { bg = it->second; bl = it->first; break; }
        } else {
            // оракул: точное попадание — где сосед умирает раньше всех позже нас (или раньше нас); иначе best fit
            int pid = (int)pos.size();
            ll my = ORD[pid];
            double best = 1e18;
            int seen = 0;
            for (auto it = byLen.lower_bound({L, -1}); it != byLen.end() && seen < 200; ++it) {
                if (it->second == ts) continue;
                ll s = it->second, l = it->first, r = l - L;
                if (bg >= 0 && r > bl - L + 8) break;
                seen++;
                ll dl = death(nb(s - 1)), dr = death(nb(s + l));
                ll soon = std::min(dl < 0 ? (1LL << 40) : dl, dr < 0 ? (1LL << 40) : dr);
                // меньше остаток — лучше; при равных — сосед умирает к нашей смерти
                double c = r * 1000.0 + (double)std::llabs(std::min<ll>(soon, 1LL << 30) - std::min<ll>(my, 1LL << 30)) / 1e6;
                if (c < best) { best = c; bg = s; bl = l; }
            }
        }
        N = 0;
        if (bg < 0) {
            ll tl = T - ts;
            if (tl < L) { N = (L - tl + 19) / 20; ll o = T; T += 20 * N; owner.resize(T, 0); freeR(o, T - o); }
            bg = tailStart(); D = bg;
        } else {
            D = bg;
            if (mode == 2 && bl > L) {   // остаток — к соседу, который умрёт раньше
                ll dl = death(nb(bg - 1)), dr = death(nb(bg + bl));
                if (dl < 0) dl = 1LL << 40; if (dr < 0) dr = 1LL << 40;
                if (dr < dl) D = bg;     // остаток справа — рядом с правым, умирающим раньше
                else D = bg + bl - L;
            }
        }
        take(bg, D, L);
        for (int i = 0; i < L; i++) owner[D + i] = (int)pos.size();
        pos.push_back(D); len.push_back(L);
    }
    ll remove(int pid) { ll d = pos[pid]; freeR(d, len[pid]); return d; }
};
