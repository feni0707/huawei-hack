// Внутрипроцессный стенд: гоняет решение (solutions/portfolio.cpp, включённый как заголовок) по тестам
// без интерактора. Считает KPI = live/trench, проверяет корректность укладок (как интерактор), меряет CPU.
//
//   sim [--pol K | --pol all | --port] [--check] [-v] [--tl] <dir|file.in>...
//     --pol K    чистая политика с номером K (см. POLS в portfolio.cpp); --pol all — матрица всех политик + oracle
//     --port     портфель (Engine с выбором политики)
//   Каждая папка — отдельный набор: строка со средним KPI (для open, где веса ≠ 1, — Σ w·KPI) и средним UB.
#define SIM
#include "../solutions/portfolio.cpp"

#include <atomic>
#include <chrono>
#include <thread>
#include <memory>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

namespace fs = std::filesystem;

struct Test {
    std::string name;
    double w = 1;
    std::vector<int> cmds;
    double ub = 0;
};

static Test loadTest(const fs::path& p) {
    Test t;
    t.name = p.filename().string();
    std::ifstream f(p);
    int n;
    f >> t.w >> n;
    t.cmds.resize(n);
    for (int i = 0; i < n; i++) f >> t.cmds[i];
    // UB: траншея >= 20*ceil(пик живых/20); KPI <= живые в конце / это
    std::vector<int> ln(1, 0);
    long long cur = 0, peak = 0;
    for (int c : t.cmds) {
        if (c > 0) { ln.push_back(c); cur += c; peak = std::max(peak, cur); }
        else { cur -= ln[-c]; }
    }
    t.ub = (double)cur / (20.0 * ((peak + 19) / 20));
    return t;
}

struct Result {
    double kpi = 0;
    long long trench = 0, live = 0;
    bool ok = true;
    std::string err;
    double ms = 0;
};

// Проверка корректности по клеткам, как в интеракторе
struct Checker {
    std::vector<int> cells;
    std::vector<int> pos, len;
    Checker() { pos.assign(1, -1); len.assign(1, 0); }
    bool place(int L, int N, int D, std::string& err) {
        for (int i = 0; i < N; i++) cells.insert(cells.end(), 20, -1);
        if (D < 0 || D + L > (int)cells.size()) { err = "out of trench"; return false; }
        int pid = (int)pos.size();
        for (int i = D; i < D + L; i++) {
            if (cells[i] != -1) { err = "collision with " + std::to_string(cells[i]); return false; }
            cells[i] = pid;
        }
        pos.push_back(D); len.push_back(L);
        return true;
    }
    bool remove(int pid, int D, std::string& err) {
        if (pos[pid] != D) { err = "wrong position of " + std::to_string(pid); return false; }
        for (int i = D; i < D + len[pid]; i++) cells[i] = -1;
        return true;
    }
};

static bool CHECK = false;
static bool DUMP = false;
static bool TRACE = false;
static bool AGG = false, LINE = false;
static bool ORACLE = false;
static std::string CSV;
static int THREADS = 0;
static void dumpHoles(const Alloc& a) {
    long long cnt[41] = {0}, cells[41] = {0};
    a.forEachHole([&](int, int len, int) { int h = std::min(len, 40); cnt[h]++; cells[h] += len; });
    long long tot = 0;
    printf("   holes (size: count/cells):");
    for (int h = 1; h <= 40; h++) if (cnt[h]) { printf(" %d%s:%lld/%lld", h, h == 40 ? "+" : "", cnt[h], cells[h]); tot += cells[h]; }
    printf("\n   total hole cells=%lld, tail=%d, trench=%d, live=%lld\n", tot, a.trench - a.tailStart, a.trench, a.live);
}

template <class F>
static Result run(const Test& t, F makeEngine) {
    Result r;
    auto t0 = std::chrono::steady_clock::now();
    auto eng = makeEngine();
    std::vector<int> death;
    if (ORACLE) {   // оракул времени смерти: предсказатель заменяется истинными данными
        int np = 0;
        for (int c : t.cmds) if (c > 0) np++;
        death.assign(np + 2, 1 << 28);
        for (size_t i = 0; i < t.cmds.size(); i++) if (t.cmds[i] < 0) death[-t.cmds[i]] = (int)i;
        oracleDeath() = &death;
    }
    Checker ck;
    size_t ci = 0;
    for (int c : t.cmds) {
        ci++;
        if (TRACE && ci % (t.cmds.size() / 20) == 0) {
            long long hc = 0, small = 0, nh = 0;
            eng->real.forEachHole([&](int, int len, int) { hc += len; nh++; if (len <= 7) small += len; });
            printf("   t=%zu trench=%d live=%lld holes=%lld (%.2f%%) nholes=%lld small<=7:%lld tail=%d\n", ci, eng->real.trench, eng->real.live, hc, 100.0 * hc / eng->real.trench, nh, small, eng->real.trench - eng->real.tailStart);
        }
        if (c > 0) {
            int N, D;
            eng->place(c, N, D);
            if (CHECK && !ck.place(c, N, D, r.err)) { r.ok = false; break; }
        } else {
            int D = eng->remove(-c);
            if (CHECK && !ck.remove(-c, D, r.err)) { r.ok = false; break; }
        }
    }
    oracleDeath() = nullptr;
    auto t1 = std::chrono::steady_clock::now();
    r.ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    r.trench = eng->real.trench;
    r.live = eng->real.live;
    r.kpi = r.ok && r.trench ? (double)r.live / r.trench : 0;
    if (DUMP) dumpHoles(eng->real);
    return r;
}

struct SetData {
    std::string name, label;
    std::vector<Test> tests;
    bool weighted = false;
};

int main(int argc, char** argv) {
    leanEnabled() = false;   // CPU-время процесса в стенде общее для всех тестов — аварийный режим выключаем (--lean включает)
    std::string mode = "port";
    int polIdx = 0;
    bool verbose = false;
    std::vector<std::string> paths;
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--pol") { std::string v = argv[++i]; if (v == "all") mode = "all"; else {
            mode = "pol"; polIdx = -1;
            for (int kk = 0; kk < NPOL; kk++) if (v == POLS[kk].name) polIdx = kk;
            if (polIdx < 0) polIdx = atoi(v.c_str());
        } }
        else if (a == "--port") mode = "port";
        else if (a == "--check") CHECK = true;
        else if (a == "-v") verbose = true;
        else if (a == "--dump") DUMP = true;
        else if (a == "--trace") TRACE = true;
        else if (a == "--agg") AGG = true;
        else if (a == "--line") { AGG = true; LINE = true; }
        else if (a == "-t") THREADS = atoi(argv[++i]);
        else if (a == "--eng") {  // --eng hyst=0.02,absMargin=20,dwell=500
            std::string v = argv[++i];
            size_t st = 0;
            while (st < v.size()) {
                size_t e = v.find(',', st);
                if (e == std::string::npos) e = v.size();
                std::string kv = v.substr(st, e - st);
                size_t q2 = kv.find('=');
                std::string f = kv.substr(0, q2);
                double val = atof(kv.substr(q2 + 1).c_str());
                if (f == "hyst") engineParams().hyst = val; else if (f == "absMargin") engineParams().absMargin = (int)val;
                else if (f == "dwell") engineParams().dwell = (int)val; else if (f == "every") engineParams().every = (int)val; else if (f == "gamma") engineParams().gamma = val; else if (f == "qwin") engineParams().qwin = (int)val; else if (f == "qq") engineParams().qq = val;
                else if (f == "hystA") engineParams().hystArena = val; else if (f == "fast") engineParams().fast = val;
                else if (f == "dwellFast") engineParams().dwellFast = (int)val; else if (f == "startCmds") engineParams().startCmds = (int)val;
                else fprintf(stderr, "unknown eng field %s\n", f.c_str());
                st = e + 1;
            }
        }
        else if (a == "--pool") {  // файл: name sel tieHigh right tailFirst js Cdiv Cmin Cmax Tsmall tuU Cmode th
            std::ifstream pf(argv[++i]);
            static std::vector<std::string> names;
            names.reserve(8192);
            POLS.clear();
            std::string line;
            while (std::getline(pf, line)) {   // 12 чисел после имени + необязательный tau*100
                std::istringstream is(line);
                std::string nm;
                int v[18] = {0};
                if (!(is >> nm)) continue;
                int cnt = 0;
                while (cnt < 18 && (is >> v[cnt])) cnt++;
                if (cnt < 12) continue;
                names.push_back(nm);
                Pol pp = {names.back().c_str(), (Sel)v[0], (bool)v[1], (bool)v[2], (bool)v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11]};
                if (cnt >= 13) pp.tau = (float)v[12] / 100.0f;
                if (cnt >= 14) pp.sw0 = (float)v[13] / 100.0f;
                if (cnt >= 15) pp.swsmall = (float)v[14] / 100.0f;
                if (cnt >= 16) pp.swsize = (float)v[15] / 100.0f;
                if (cnt >= 17) pp.wkap = (float)v[16] / 100.0f;
                POLS.push_back(pp);
            }
        }
        else if (a == "--csv") CSV = argv[++i];
        else if (a == "--dbgswitch") dbgSwitch() = true;
        else if (a == "--oracle") ORACLE = true;
        else if (a == "--lean") leanEnabled() = true;
        else if (a == "--fam") g_famMask = (int)strtol(argv[++i], nullptr, 0);
        else if (a == "--nh") g_nhLimit = atoi(argv[++i]);
        else if (a == "--pset") {  // --pset idx:field=val[,field=val]  (js, Cdiv, Cmin, Cmax, Tsmall, right, tieHigh, tailFirst, sel)
            std::string v = argv[++i];
            size_t c = v.find(':');
            std::string nm = v.substr(0, c);
            int idx = -1;
            for (int kk = 0; kk < NPOL; kk++) if (nm == POLS[kk].name) idx = kk;
            if (idx < 0) idx = atoi(nm.c_str());
            std::string rest = v.substr(c + 1);
            size_t st = 0;
            while (st < rest.size()) {
                size_t e = rest.find(',', st);
                if (e == std::string::npos) e = rest.size();
                std::string kv = rest.substr(st, e - st);
                size_t q2 = kv.find('=');
                std::string f = kv.substr(0, q2);
                int val = atoi(kv.substr(q2 + 1).c_str());
                Pol& pp = POLS[idx];
                if (f == "js") pp.js = val; else if (f == "Cdiv") pp.Cdiv = val; else if (f == "Cmin") pp.Cmin = val;
                else if (f == "Cmax") pp.Cmax = val; else if (f == "Tsmall") pp.Tsmall = val; else if (f == "tuU") pp.tuU = val; else if (f == "Cmode") pp.Cmode = val; else if (f == "th") pp.th = val; else if (f == "right") pp.right = val;
                else if (f == "tieHigh") pp.tieHigh = val; else if (f == "tailFirst") pp.tailFirst = val; else if (f == "sel") pp.sel = (Sel)val;
                else if (f == "tau") pp.tau = (float)val / 100.0f; else if (f == "sw0") pp.sw0 = (float)val / 100.0f;
                else if (f == "swsmall") pp.swsmall = (float)val / 100.0f; else if (f == "swsize") pp.swsize = (float)val / 100.0f;
                else if (f == "wkap") pp.wkap = (float)val / 100.0f;
                else fprintf(stderr, "unknown field %s\n", f.c_str());
                st = e + 1;
            }
        }
        else if (a == "--dump") DUMP = true;
        else paths.push_back(a);
    }
    std::vector<SetData> sets;
    for (auto& p : paths) {
        SetData sd;
        sd.name = fs::path(p).parent_path().filename().string() + "/" + fs::path(p).filename().string();
        if (fs::is_directory(p)) {
            sd.name = fs::path(p).filename().string();
            if (sd.name.empty()) sd.name = fs::path(p).parent_path().filename().string();
            sd.label = fs::path(p).parent_path().filename().string() + "/" + sd.name;
            std::vector<fs::path> files;
            for (auto& e : fs::directory_iterator(p)) if (e.path().extension() == ".in") files.push_back(e.path());
            std::sort(files.begin(), files.end(), [](const fs::path& a, const fs::path& b) {
                std::string x = a.stem().string(), y = b.stem().string();
                if (x.size() != y.size()) return x.size() < y.size();
                return x < y;
            });
            for (auto& f : files) sd.tests.push_back(loadTest(f));
        } else {
            sd.tests.push_back(loadTest(p));
        }
        for (auto& t : sd.tests) if (t.w != 1) sd.weighted = true;
        sets.push_back(std::move(sd));
    }

    auto agg = [&](const SetData& sd, const std::vector<double>& k) {
        double s = 0, sw = 0;
        for (size_t i = 0; i < sd.tests.size(); i++) {
            double w = sd.weighted ? sd.tests[i].w : 1;
            s += k[i] * w; sw += w;
        }
        return sd.weighted ? s : s / std::max(1.0, sw);
    };

    if (AGG) {
        std::vector<SetData> merged;
        for (auto& sd : sets) {
            bool found = false;
            for (auto& m : merged) if (m.name == sd.name) { for (auto& t : sd.tests) m.tests.push_back(t); found = true; break; }
            if (!found) merged.push_back(sd);
        }
        sets = std::move(merged);
    }
    if (mode == "all") {
        // матрица: каждая чистая политика на каждом тесте (параллельно), oracle = лучшая на каждом тесте, wins = сколько раз политика лучшая
        int NP = NPOL;
        struct J { size_t set, idx; int k; double kpi; };
        std::vector<J> jobs;
        for (size_t si = 0; si < sets.size(); si++) for (size_t i = 0; i < sets[si].tests.size(); i++) for (int k = 0; k < NP; k++) jobs.push_back({si, i, k, 0});
        std::atomic<size_t> next(0);
        int nth = (int)std::max(1u, std::thread::hardware_concurrency());
        std::vector<std::thread> pool;
        for (int ti = 0; ti < nth; ti++) pool.emplace_back([&] {
            for (;;) {
                size_t j = next++;
                if (j >= jobs.size()) break;
                Result r = run(sets[jobs[j].set].tests[jobs[j].idx], [&] { return std::unique_ptr<Engine>(new Engine(jobs[j].k)); });
                jobs[j].kpi = r.kpi;
                if (!r.ok) fprintf(stderr, "INVALID pol %d%s", jobs[j].k, "\n");
            }
        });
        for (auto& th : pool) th.join();
        if (!CSV.empty()) {
            FILE* cf = fopen(CSV.c_str(), "w");
            fprintf(cf, "set,test,w,pol,kpi%s", "\n");
            for (auto& jb : jobs) fprintf(cf, "%s,%s,%g,%s,%.6f%s", sets[jb.set].label.c_str(), sets[jb.set].tests[jb.idx].name.c_str(), sets[jb.set].tests[jb.idx].w, POLS[jb.k].name, jb.kpi, "\n");
            fclose(cf);
            return 0;
        }
        printf("%-16s", "set");
        for (int k = 0; k < NP; k++) printf(" %8s", POLS[k].name);
        printf(" %8s %8s\n", "oracle", "UB");
        std::vector<long long> wins(NP, 0);
        long long ntests = 0;
        for (size_t si = 0; si < sets.size(); si++) {
            size_t nt = sets[si].tests.size();
            std::vector<std::vector<double>> K(NP, std::vector<double>(nt, 0));
            std::vector<double> orc(nt, 0), ub(nt);
            for (auto& jb : jobs) if (jb.set == si) { K[jb.k][jb.idx] = jb.kpi; orc[jb.idx] = std::max(orc[jb.idx], jb.kpi); }
            for (size_t i = 0; i < nt; i++) {
                ub[i] = sets[si].tests[i].ub;
                ntests++;
                for (int k = 0; k < NP; k++) if (K[k][i] >= orc[i] - 1e-9) wins[k]++;
            }
            printf("%-16s", sets[si].name.c_str());
            for (int k = 0; k < NP; k++) printf(" %8.4f", agg(sets[si], K[k]));
            printf(" %8.4f %8.4f\n", agg(sets[si], orc), agg(sets[si], ub));
        }
        printf("%-16s", "wins/tests");
        for (int k = 0; k < NP; k++) printf(" %8lld", wins[k]);
        printf("  of %lld\n", ntests);
        return 0;
    }

    // параллельный прогон всех тестов
    struct Job { size_t set, idx; Result r; };
    std::vector<Job> jobs;
    for (size_t si = 0; si < sets.size(); si++) for (size_t i = 0; i < sets[si].tests.size(); i++) jobs.push_back({si, i, Result()});
    std::atomic<size_t> next(0);
    int nth = (int)std::max(1u, std::thread::hardware_concurrency());
    if (THREADS > 0) nth = THREADS;
    std::vector<std::thread> pool;
    for (int ti = 0; ti < nth; ti++) pool.emplace_back([&] {
        for (;;) {
            size_t j = next++;
            if (j >= jobs.size()) break;
            const Test& t = sets[jobs[j].set].tests[jobs[j].idx];
            jobs[j].r = mode == "pol" ? run(t, [&] { return std::unique_ptr<Engine>(new Engine(polIdx)); })
                                      : run(t, [&] { return std::unique_ptr<Engine>(new Engine(-1)); });
        }
    });
    for (auto& th : pool) th.join();
    std::vector<std::vector<double>> K(sets.size());
    std::vector<double> maxms(sets.size(), 0);
    for (size_t si = 0; si < sets.size(); si++) K[si].assign(sets[si].tests.size(), 0);
    for (auto& jb : jobs) {
        const Test& t = sets[jb.set].tests[jb.idx];
        if (!jb.r.ok) fprintf(stderr, "INVALID %s test %s: %s%s", sets[jb.set].name.c_str(), t.name.c_str(), jb.r.err.c_str(), "\n");
        K[jb.set][jb.idx] = jb.r.kpi;
        maxms[jb.set] = std::max(maxms[jb.set], jb.r.ms);
        if (verbose) printf("  %-16s %-8s w=%-4g kpi=%.4f ub=%.4f trench=%lld live=%lld  %.1fms\n", sets[jb.set].name.c_str(), t.name.c_str(), t.w, jb.r.kpi, t.ub, jb.r.trench, jb.r.live, jb.r.ms);
    }
    if (LINE) {
        for (size_t si = 0; si < sets.size(); si++) printf("%s=%.4f ", sets[si].name.c_str(), agg(sets[si], K[si]));
        printf("\n");
        return 0;
    }
    for (size_t si = 0; si < sets.size(); si++) {
        std::vector<double> ub(sets[si].tests.size());
        for (size_t i = 0; i < ub.size(); i++) ub[i] = sets[si].tests[i].ub;
        printf("%-16s %s=%.4f  UB=%.4f  maxms=%.0f\n", sets[si].name.c_str(), sets[si].weighted ? "sum" : "avg", agg(sets[si], K[si]), agg(sets[si], ub), maxms[si]);
    }
    return 0;
}
