// Best-fit allocator with coalescing free list.
// Trench grows only when no free gap fits; growth is minimal and reuses the tail gap.
#include <bits/stdc++.h>
using namespace std;

int main() {
    long long trenchLen = 0;
    map<long long, long long> freeByStart;          // start -> length
    set<pair<long long, long long>> freeBySize;     // (length, start)
    vector<long long> posOf(1, -1), lenOf(1, 0);    // by PID

    auto addFree = [&](long long s, long long l) {
        if (l <= 0) return;
        auto it = freeByStart.lower_bound(s);
        if (it != freeByStart.end() && it->first == s + l) {
            freeBySize.erase({it->second, it->first});
            l += it->second;
            it = freeByStart.erase(it);
        }
        if (it != freeByStart.begin()) {
            auto pv = prev(it);
            if (pv->first + pv->second == s) {
                freeBySize.erase({pv->second, pv->first});
                s = pv->first;
                l += pv->second;
                freeByStart.erase(pv);
            }
        }
        freeByStart[s] = l;
        freeBySize.insert({l, s});
    };
    auto removeFree = [&](long long s, long long l) {
        freeBySize.erase({l, s});
        freeByStart.erase(s);
    };

    long long x;
    while (scanf("%lld", &x) == 1) {
        if (x == 0) break;
        if (x > 0) {
            long long L = x, N = 0, D;
            auto it = freeBySize.lower_bound({L, -1});
            if (it != freeBySize.end()) {
                long long gl = it->first, gs = it->second;
                D = gs;
                removeFree(gs, gl);
                addFree(gs + L, gl - L);
            } else {
                long long tailS = trenchLen, tailL = 0;
                if (!freeByStart.empty()) {
                    auto last = prev(freeByStart.end());
                    if (last->first + last->second == trenchLen) {
                        tailS = last->first;
                        tailL = last->second;
                    }
                }
                N = (L - tailL + 19) / 20;
                if (tailL > 0) removeFree(tailS, tailL);
                trenchLen += N * 20;
                D = tailS;
                addFree(D + L, trenchLen - (D + L));
            }
            posOf.push_back(D);
            lenOf.push_back(L);
            printf("%lld %lld\n", N, D);
        } else {
            long long pid = -x;
            long long s = posOf[pid], l = lenOf[pid];
            printf("%lld\n", s);
            addFree(s, l);
        }
        fflush(stdout);
    }
    return 0;
}
