// Внутрипроцессный стенд: гоняет решение (solutions/portfolio.cpp, включённый как заголовок) по тестам
// без интерактора. Считает KPI = live/trench, проверяет корректность укладок (как интерактор), меряет CPU.
//
//   sim [--pol K | --pol all | --port] [--check] [-v] [--tl] <dir|file.in>...
//     --pol K    чистая политика с номером K (см. POLS в portfolio.cpp); --pol all — матрица всех политик + oracle
//     --port     портфель (Engine с выбором политики)
//   Каждая папка — отдельный набор: строка со средним KPI (для open, где веса ≠ 1, — Σ w·KPI) и средним UB.
#define SIM
#include "../solutions/portfolio.cpp"

#include <chrono>
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
static void dumpHoles(const Alloc& a) {
    long long cnt[41] = {0}, cells[41] = {0};
    for (auto& kv : a.byAddr) { int h = std::min(kv.second.len, 40); cnt[h]++; cells[h] += kv.second.len; }
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
    Checker ck;
    size_t ci = 0;
    for (int c : t.cmds) {
        ci++;
        if (TRACE && ci % (t.cmds.size() / 20) == 0) {
            long long hc = 0, small = 0, nh = 0;
            for (auto& kv : eng->real.byAddr) { hc += kv.second.len; nh++; if (kv.second.len <= 7) small += kv.second.len; }
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
    auto t1 = std::chrono::steady_clock::now();
    r.ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    r.trench = eng->real.trench;
    r.live = eng->real.live;
    r.kpi = r.ok && r.trench ? (double)r.live / r.trench : 0;
    if (DUMP) dumpHoles(eng->real);
    return r;
}

struct SetData {
    std::string name;
    std::vector<Test> tests;
    bool weighted = false;
};

int main(int argc, char** argv) {
    std::string mode = "port";
    int polIdx = 0;
    bool verbose = false;
    std::vector<std::string> paths;
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--pol") { std::string v = argv[++i]; if (v == "all") mode = "all"; else { mode = "pol"; polIdx = atoi(v.c_str()); } }
        else if (a == "--port") mode = "port";
        else if (a == "--check") CHECK = true;
        else if (a == "-v") verbose = true;
        else if (a == "--dump") DUMP = true;
        else if (a == "--trace") TRACE = true;
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

    if (mode == "all") {
        printf("%-16s", "set");
        for (int k = 0; k < NPOL; k++) printf(" %8s", POLS[k].name);
        printf(" %8s %8s\n", "oracle", "UB");
        std::vector<double> colTot(NPOL + 2, 0);
        for (auto& sd : sets) {
            std::vector<std::vector<double>> K(NPOL, std::vector<double>(sd.tests.size()));
            std::vector<double> orc(sd.tests.size(), 0), ub(sd.tests.size());
            for (int k = 0; k < NPOL; k++)
                for (size_t i = 0; i < sd.tests.size(); i++) {
                    Result r = run(sd.tests[i], [&] { return std::unique_ptr<Engine>(new Engine(k)); });
                    K[k][i] = r.kpi;
                    orc[i] = std::max(orc[i], r.kpi);
                    if (!r.ok) fprintf(stderr, "INVALID %s pol %s test %s: %s\n", sd.name.c_str(), POLS[k].name, sd.tests[i].name.c_str(), r.err.c_str());
                }
            for (size_t i = 0; i < sd.tests.size(); i++) ub[i] = sd.tests[i].ub;
            printf("%-16s", sd.name.c_str());
            for (int k = 0; k < NPOL; k++) printf(" %8.4f", agg(sd, K[k]));
            printf(" %8.4f %8.4f\n", agg(sd, orc), agg(sd, ub));
        }
        return 0;
    }

    double grand = 0;
    for (auto& sd : sets) {
        std::vector<double> k(sd.tests.size()), ub(sd.tests.size());
        double maxms = 0;
        for (size_t i = 0; i < sd.tests.size(); i++) {
            Result r = mode == "pol" ? run(sd.tests[i], [&] { return std::unique_ptr<Engine>(new Engine(polIdx)); })
                                     : run(sd.tests[i], [&] { return std::unique_ptr<Engine>(new Engine(-1)); });
            if (!r.ok) fprintf(stderr, "INVALID %s test %s: %s\n", sd.name.c_str(), sd.tests[i].name.c_str(), r.err.c_str());
            k[i] = r.kpi; ub[i] = sd.tests[i].ub;
            maxms = std::max(maxms, r.ms);
            if (verbose) printf("  %-16s %-8s w=%-4g kpi=%.4f ub=%.4f trench=%lld live=%lld  %.1fms\n", sd.name.c_str(),
                                sd.tests[i].name.c_str(), sd.tests[i].w, r.kpi, sd.tests[i].ub, r.trench, r.live, r.ms);
        }
        printf("%-16s %s=%.4f  UB=%.4f  maxms=%.0f\n", sd.name.c_str(), sd.weighted ? "sum" : "avg", agg(sd, k), agg(sd, ub), maxms);
    }
    return 0;
}
