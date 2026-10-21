// Оффлайн-оракул: оценка того, чего можно достичь, ЗНАЯ будущее (потолок для предсказательных стратегий).
//
// Раскладка «обратного времени»: идём по командам с конца. Трубы, пережившие весь тест, укладываются
// сразу и плотно (порядок — по PID); далее обратная по времени команда «уложить L» освобождает трубу,
// а «удалить K» — укладывает её (политика best fit из solutions/variant_a.cpp). Так получается допустимая
// раскладка для известной заранее последовательности; траншея = максимальная протяжённость.
// KPI оракула = (живые в конце) / (20·ceil(макс. протяжённость / 20)).
//
// Сборка: c++ -std=c++17 -O2 -o build/oracle scripts/oracle.cpp
// Запуск: build/oracle test1.in test2.in ...   → строка на тест: имя  kpi_forward  kpi_oracle  ub
#define main variant_a_main
#include "../solutions/variant_a.cpp"
#undef main

#include <cstdio>
#include <vector>
using namespace std;

static double runForward(const vector<ll>& cmds, ll& liveFinal, ll& peak) {
    Params P;
    Allocator A(P);
    for (ll c : cmds) {
        if (c > 0) A.place(c);
        else A.remove(-c);
    }
    liveFinal = A.live;
    peak = A.peak;
    return A.trench() ? (double)A.live / A.trench() : 1.0;
}

static double runOracle(const vector<ll>& cmds, int order) {
    int n = cmds.size();
    // длины и «выжившие»
    vector<ll> len(1, 0);
    vector<char> removed(1, 0);
    for (ll c : cmds) {
        if (c > 0) { len.push_back(c); removed.push_back(0); }
        else removed[-c] = 1;
    }
    int m = len.size() - 1;
    Params P;
    Allocator A(P);
    vector<ll> bpid(m + 1, 0);   // PID в обратной раскладке
    ll live = 0;
    vector<int> surv;
    for (int p = 1; p <= m; p++)
        if (!removed[p]) surv.push_back(p);
    if (order == 1) {            // выжившие: сначала длинные
        stable_sort(surv.begin(), surv.end(), [&](int a, int b) { return len[a] > len[b]; });
    }
    ll cnt = 0;
    for (int p : surv) {
        A.place(len[p]);
        bpid[p] = ++cnt;
        live += len[p];
    }
    ll T = A.trench();
    int pid = m;
    for (int i = n - 1; i >= 0; i--) {
        ll c = cmds[i];
        if (c > 0) {
            if (bpid[pid]) A.remove(bpid[pid]);
            pid--;
        } else {
            ll p = -c;
            A.place(len[p]);
            bpid[p] = ++cnt;
            if (A.trench() > T) T = A.trench();
        }
    }
    return T ? (double)live / T : 1.0;
}

int main(int argc, char** argv) {
    for (int a = 1; a < argc; a++) {
        FILE* f = fopen(argv[a], "r");
        if (!f) return 2;
        double w;
        ll n;
        if (fscanf(f, "%lf %lld", &w, &n) != 2) return 2;
        vector<ll> cmds(n);
        for (ll i = 0; i < n; i++)
            if (fscanf(f, "%lld", &cmds[i]) != 1) return 2;
        fclose(f);
        ll live, peak;
        double kf = runForward(cmds, live, peak);
        double ko = max(runOracle(cmds, 0), runOracle(cmds, 1));
        ll ubT = 20 * ((peak + 19) / 20);
        printf("%s %.6f %.6f %.6f\n", argv[a], kf, ko, ubT ? (double)live / ubT : 1.0);
    }
    return 0;
}
