// LNS: destroy (3 versions) + repair (weighted 2-regret, greedy cycle), with and without local search
// Build: g++ -O2 -std=c++17 -o lns lns.cpp
// Usage: ./lns [--data-dir DIR] [--out DIR] [--runs N] [--seed S] [--frac F] [--msls FILE]
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
typedef vector<int> Tour;

// ---------------------------------------------------------------- instance
struct Instance {
    string name;
    int n = 0, k = 0;
    vector<double> x, y;
    vector<int> cost;
    vector<vector<int>> D;
    vector<int> F;   // F[a * n + b]
};

Instance load_instance(const string &dir, const string &name) {
    Instance in;
    in.name = name;
    ifstream f(dir + "/" + name + ".csv");
    if (!f) { cerr << "cannot open " << dir << "/" << name << ".csv\n"; exit(1); }
    string line;
    while (getline(f, line)) {
        if (line.empty()) continue;
        for (char &c : line) if (c == ';') c = ' ';
        istringstream ss(line);
        double a, b; int c;
        ss >> a >> b >> c;
        in.x.push_back(a); in.y.push_back(b); in.cost.push_back(c);
    }
    in.n = in.x.size();
    in.k = (in.n + 1) / 2;
    in.D.assign(in.n, vector<int>(in.n, 0));
    for (int i = 0; i < in.n; i++)
        for (int j = 0; j < in.n; j++)
            in.D[i][j] = (int)floor(sqrt((in.x[i] - in.x[j]) * (in.x[i] - in.x[j]) +
                                         (in.y[i] - in.y[j]) * (in.y[i] - in.y[j])) + 0.5);
    in.F.resize(in.n * in.n);
    for (int i = 0; i < in.n; i++)
        for (int j = 0; j < in.n; j++) in.F[i * in.n + j] = in.D[i][j];
    return in;
}

long long objective(const Instance &in, const Tour &t) {
    long long s = 0;
    int k = t.size();
    for (int i = 0; i < k; i++) s += in.D[t[i]][t[(i + 1) % k]] + in.cost[t[i]];
    return s;
}

// independent checker: uses coordinates, not D
long long check(const Instance &in, const Tour &t) {
    vector<char> seen(in.n, 0);
    if ((int)t.size() != in.k) { cerr << "checker: wrong size\n"; exit(1); }
    long long s = 0;
    int k = t.size();
    for (int i = 0; i < k; i++) {
        if (seen[t[i]]) { cerr << "checker: repeated node\n"; exit(1); }
        seen[t[i]] = 1;
        int a = t[i], b = t[(i + 1) % k];
        double dx = in.x[a] - in.x[b], dy = in.y[a] - in.y[b];
        s += (int)(sqrt(dx * dx + dy * dy) + 0.5) + in.cost[a];
    }
    return s;
}

typedef mt19937_64 Rng;
inline int rnd(Rng &r, int m) { return (int)(((r() >> 32) * (unsigned long long)m) >> 32); }

Tour random_solution(const Instance &in, Rng &r) {
    Tour p(in.n);
    for (int i = 0; i < in.n; i++) p[i] = i;
    for (int i = 0; i < in.k; i++) swap(p[i], p[i + rnd(r, in.n - i)]);
    p.resize(in.k);
    return p;
}

// ---------------------------------------------------------------- steepest with LM
// type 0: node exchange (v, u, p, nx), type 1: edges exchange (a -> an, b -> bn)
struct LMSteepest {
    struct Move { int d, type, a, an, b, bn; };
    const Instance &in;
    int K, n;
    const int *F, *C;
    Tour t;
    vector<int> pos;   // position in t, -1 if not selected
    // heap of keys (delta << 32 | index in moves)
    vector<Move> moves;
    vector<long long> LM;
    vector<long long> later;   // SKIP moves
    bool heap = false;
    long long evals = 0, applied = 0;

    LMSteepest(const Instance &in_) : in(in_) { K = in.k; n = in.n; F = in.F.data(); C = in.cost.data(); }

    inline int nxt(int i) const { return i + 1 < K ? i + 1 : 0; }
    inline int prv(int i) const { return i ? i - 1 : K - 1; }
    inline int succ(int v) const { return t[nxt(pos[v])]; }
    inline int pred(int v) const { return t[prv(pos[v])]; }

    // edges exchange: remove (a -> an), (b -> bn), add (a, b), (an, bn)
    inline void eval_edges(int a, int an, int b, int bn) {
        int d = F[a * n + b] + F[an * n + bn] - F[a * n + an] - F[b * n + bn];
        evals++;
        if (d < 0) push({d, 1, a, an, b, bn});
    }
    // both versions for edges starting at i and j
    inline void eval_edge_pair(int i, int j) {
        int a = t[i], an = t[nxt(i)], b = t[j], bn = t[nxt(j)];
        eval_edges(a, an, b, bn);
        eval_edges(a, an, bn, b);
    }
    // edge starting at i x all other edges
    void eval_edge_all(int i, int skip = -1) {
        int a = t[i], an = t[nxt(i)], ip = prv(i), in_ = nxt(i);
        const int *Fa = F + a * n, *Fan = F + an * n;
        int d_a_an = Fa[an];
        for (int j = 0; j < K; j++) {
            if (j == i || j == in_ || j == ip || j == skip) continue;
            int b = t[j], bn = t[j + 1 < K ? j + 1 : 0];
            int base = -d_a_an - F[b * n + bn];
            int d1 = base + Fa[b] + Fan[bn];   // remove (a -> an), (b -> bn), add (a, b), (an, bn)
            int d2 = base + Fa[bn] + Fan[b];   // second edge reversed
            evals += 2;
            if (d1 < 0) push({d1, 1, a, an, b, bn});
            if (d2 < 0) push({d2, 1, a, an, bn, b});
        }
    }
    // node exchange: v (selected) replaced by u (not selected)
    inline void eval_node(int v, int u) {
        int p = pred(v), nx = succ(v);
        int d = F[p * n + u] + F[u * n + nx] - F[p * n + v] - F[v * n + nx] + C[u] - C[v];
        evals++;
        if (d < 0) push({d, 0, v, u, p, nx});
    }
    void eval_node_all(int v) {
        int p = pred(v), nx = succ(v);
        const int *Fp = F + p * n, *Fn = F + nx * n;
        int base = -Fp[v] - Fn[v] - C[v];
        for (int u = 0; u < n; u++) {
            if (pos[u] >= 0) continue;
            int d = base + Fp[u] + Fn[u] + C[u];
            evals++;
            if (d < 0) push({d, 0, v, u, p, nx});
        }
    }
    // min-heap
    static bool worse(long long x, long long y) { return y < x; }
    inline void push_key(long long key) {
        LM.push_back(key);
        if (heap) push_heap(LM.begin(), LM.end(), worse);
    }
    inline void push(const Move &m) {
        push_key((long long)m.d * (1LL << 32) + (long long)moves.size());
        moves.push_back(m);
    }
    inline long long pop() {
        pop_heap(LM.begin(), LM.end(), worse);
        long long key = LM.back(); LM.pop_back();
        return key;
    }
    // +1: a -> b is an edge, -1: b -> a is an edge, 0: no edge
    inline int dir(int a, int b) const {
        if (pos[a] < 0 || pos[b] < 0) return 0;
        if (succ(a) == b) return 1;
        if (pred(a) == b) return -1;
        return 0;
    }

    void apply_edges(int i, int j) {   // remove edges starting at positions i and j
        if (i > j) swap(i, j);
        reverse(t.begin() + i + 1, t.begin() + j + 1);
        for (int q = i + 1; q <= j; q++) pos[t[q]] = q;
        // new edges start at i and j
        eval_edge_all(i);
        eval_edge_all(j, i);
        for (int q : {i, nxt(i), j, nxt(j)}) eval_node_all(t[q]);
    }
    void apply_node(int v, int u) {
        int i = pos[v];
        t[i] = u; pos[u] = i; pos[v] = -1;
        // new edges start at i-1 and i
        eval_edge_all(prv(i));
        eval_edge_all(i, prv(i));
        for (int q : {prv(i), i, nxt(i)}) eval_node_all(t[q]);
        for (int q = 0; q < K; q++)   // v is outside now
            if (q != prv(i) && q != i && q != nxt(i)) eval_node(t[q], v);
    }

    Tour run(const Tour &start) {
        t = start; evals = 0; applied = 0;
        LM.clear(); moves.clear(); heap = false;
        pos.assign(n, -1);
        for (int i = 0; i < K; i++) pos[t[i]] = i;
        // initial solution: all moves
        for (int i = 0; i < K; i++)
            for (int j = i + 2; j < K; j++)
                if (!(i == 0 && j == K - 1)) eval_edge_pair(i, j);
        for (int i = 0; i < K; i++) eval_node_all(t[i]);
        make_heap(LM.begin(), LM.end(), worse);
        heap = true;

        while (true) {
            bool found = false;
            Move m{};
            later.clear();
            while (!LM.empty()) {
                long long key = pop();
                m = moves[key & 0xffffffffLL];
                if (m.type == 0) {
                    int v = m.a, u = m.an, p = m.b, nx = m.bn;
                    if (pos[v] >= 0 && pos[u] < 0 &&
                        ((pred(v) == p && succ(v) == nx) || (pred(v) == nx && succ(v) == p))) {
                        found = true; break;
                    }
                } else {
                    int d1 = dir(m.a, m.an), d2 = dir(m.b, m.bn);
                    if (d1 == 0 || d2 == 0) continue;   // REMOVE
                    if (d1 != d2) { later.push_back(key); continue; }   // SKIP
                    found = true; break;
                }
            }
            for (long long x : later) push_key(x);
            if (!found) break;
            applied++;
#ifdef CHECK_DELTA
            long long before = objective(in, t);
#endif
            if (m.type == 0) apply_node(m.a, m.an);
            else if (dir(m.a, m.an) == 1) apply_edges(pos[m.a], pos[m.b]);
            else apply_edges(pos[m.an], pos[m.bn]);   // both reversed
#ifdef CHECK_DELTA
            if (objective(in, t) - before != m.d) { cerr << "delta mismatch\n"; exit(1); }
#endif
        }
        return t;
    }
};

// ---------------------------------------------------------------- destroy
// all versions remove m = round(frac * K) nodes, the rest keeps its order
Tour destroy_nodes(const Tour &x, int m, Rng &r) {          // m random nodes
    Tour t = x;
    for (int q = 0; q < m; q++) {
        int i = rnd(r, t.size());
        t.erase(t.begin() + i);
    }
    return t;
}

Tour destroy_subpath(const Tour &x, int m, Rng &r) {        // one subpath of m nodes
    int K = x.size(), s = rnd(r, K);
    Tour t;
    for (int q = m; q < K; q++) t.push_back(x[(s + q) % K]);
    return t;
}

// random split of total into parts >= 1
vector<int> split(int total, int parts, Rng &r) {
    vector<int> cut;
    for (int i = 1; i < total; i++) cut.push_back(i);
    for (int i = 0; i < parts - 1; i++) swap(cut[i], cut[i + rnd(r, cut.size() - i)]);
    cut.resize(parts - 1);
    sort(cut.begin(), cut.end());
    vector<int> res;
    int prev = 0;
    for (int c : cut) { res.push_back(c - prev); prev = c; }
    res.push_back(total - prev);
    return res;
}

Tour destroy_subpaths(const Tour &x, int m, Rng &r) {       // 3-5 subpaths, m nodes together
    int K = x.size(), p = 3 + rnd(r, 3), s = rnd(r, K);
    vector<int> rem = split(m, p, r), keep = split(K - m, p, r);   // kept parts between them
    Tour t;
    int i = 0;
    for (int q = 0; q < p; q++) {
        i += rem[q];
        for (int e = 0; e < keep[q]; e++) t.push_back(x[(s + i++) % K]);
    }
    return t;
}

// ---------------------------------------------------------------- repair
// weighted 2-regret greedy cycle (best constructive heuristic from assignment 2):
// insert node with max  regret - 2 * best, regret = 2nd best insertion - best insertion
Tour repair(const Instance &in, Tour t) {
    const int n = in.n, *F = in.F.data(), *C = in.cost.data();
    // cycle as a successor array; for every free node u: 2 best insertion edges (start nodes e1, e2)
    vector<int> nx(n, -1), b1(n), b2(n), e1(n), e2(n);
    int m = t.size();
    for (int i = 0; i < m; i++) nx[t[i]] = t[i + 1 < m ? i + 1 : 0];
    vector<int> fr;
    for (int u = 0; u < n; u++) if (nx[u] < 0) fr.push_back(u);
    auto ins = [&](int u, int a) { int b = nx[a]; return F[u * n + a] + F[u * n + b] - F[a * n + b]; };
    auto consider = [&](int u, int a) {
        int d = ins(u, a);
        if (d < b1[u]) { b2[u] = b1[u]; e2[u] = e1[u]; b1[u] = d; e1[u] = a; }
        else if (d < b2[u]) { b2[u] = d; e2[u] = a; }
    };
    auto full = [&](int u) {
        b1[u] = b2[u] = INT_MAX; e1[u] = e2[u] = -1;
        int a = t[0];
        for (int i = 0; i < m; i++, a = nx[a]) consider(u, a);
        if (m < 2) { b2[u] = b1[u]; e2[u] = e1[u]; }
    };
    for (int u : fr) full(u);
    while (m < in.k) {
        long long bestScore = LLONG_MIN; int bi = -1, bbest = 0;
        for (int i = 0; i < (int)fr.size(); i++) {
            int u = fr[i], best = b1[u] + C[u];
            long long score = (long long)(b2[u] - b1[u]) - 2LL * best;
            if (score > bestScore || (score == bestScore && best < bbest)) { bestScore = score; bi = i; bbest = best; }
        }
        int w = fr[bi], a = e1[w], b = nx[a];
        fr[bi] = fr.back(); fr.pop_back();
        nx[a] = w; nx[w] = b; m++;
        // edge (a, b) replaced by (a, w), (w, b)
        for (int u : fr) {
            if (e1[u] == a || e2[u] == a) full(u);
            else { consider(u, a); consider(u, w); }
        }
    }
    Tour res;
    int v = t[0];
    for (int i = 0; i < m; i++, v = nx[v]) res.push_back(v);
    return res;
}

// ---------------------------------------------------------------- experiment
// average MSLS time per instance from assignment 6
map<string, double> msls_time(const string &path) {
    map<string, double> sum, cnt;
    ifstream f(path);
    if (!f) { cerr << "cannot open " << path << "\n"; exit(1); }
    string line;
    getline(f, line);
    while (getline(f, line)) {
        vector<string> c;
        string cur; bool q = false;
        for (char ch : line) {
            if (ch == '"') q = !q;
            else if (ch == ',' && !q) { c.push_back(cur); cur.clear(); }
            else cur += ch;
        }
        c.push_back(cur);
        if (c[1] == "MSLS") { sum[c[0]] += atof(c[4].c_str()); cnt[c[0]] += 1; }
    }
    for (auto &p : sum) p.second /= cnt[p.first];
    return sum;
}

int main(int argc, char **argv) {
    string dataDir = "../data", outDir = "results", mslsFile = "../lab6/results/raw.csv";
    int runs = 20; double frac = 0.3;
    unsigned long long seed = 42;
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        if (a == "--data-dir" && i + 1 < argc) dataDir = argv[++i];
        else if (a == "--out" && i + 1 < argc) outDir = argv[++i];
        else if (a == "--runs" && i + 1 < argc) runs = atoi(argv[++i]);
        else if (a == "--seed" && i + 1 < argc) seed = strtoull(argv[++i], 0, 10);
        else if (a == "--frac" && i + 1 < argc) frac = atof(argv[++i]);
        else if (a == "--msls" && i + 1 < argc) mslsFile = argv[++i];
        else { cerr << "usage: lns [--data-dir DIR] [--out DIR] [--runs N] [--seed S] [--frac F] [--msls FILE]\n"; return 1; }
    }
    map<string, double> limit = msls_time(mslsFile);
    string cmd = "mkdir -p '" + outDir + "'";
    if (system(cmd.c_str()) != 0) return 1;
    ofstream raw(outDir + "/raw.csv"), best(outDir + "/best.txt");
    raw << "instance,method,run,value,time_us,iters\n";
    typedef chrono::steady_clock Clock;
    const char *dnames[3] = {"random nodes", "single subpath", "3-5 subpaths"};

    for (string name : {"TSPA", "TSPB"}) {
        Instance in = load_instance(dataDir, name);
        LMSteepest ls(in);
        int m = (int)(frac * in.k + 0.5);
        double T = limit[name];
        for (int useLS = 1; useLS >= 0; useLS--)
        for (int dv = 0; dv < 3; dv++) {
            string method = string("LNS, ") + dnames[dv] + (useLS ? ", with LS" : ", without LS");
            Tour bestTour; long long bestVal = LLONG_MAX;
            double sv = 0, st = 0, si = 0; long long mn = LLONG_MAX, mx = 0;
            for (int r = 0; r < runs; r++) {
                Rng g(seed * 1000003ULL + r * 7919ULL + dv * 31 + useLS);
                auto t0 = Clock::now();
                Tour x = ls.run(random_solution(in, g));
                long long fx = objective(in, x);
                int iters = 0;
                while (chrono::duration<double, micro>(Clock::now() - t0).count() < T) {
                    Tour y = dv == 0 ? destroy_nodes(x, m, g) : dv == 1 ? destroy_subpath(x, m, g)
                                                                        : destroy_subpaths(x, m, g);
                    y = repair(in, y);
                    if (useLS) y = ls.run(y);
                    long long fy = objective(in, y);
                    iters++;
                    if (fy < fx) { fx = fy; x = y; }
                }
                double us = chrono::duration<double, micro>(Clock::now() - t0).count();
                raw << name << ",\"" << method << "\"," << r << "," << fx << "," << us << "," << iters << "\n";
                if (fx < bestVal) { bestVal = fx; bestTour = x; }
                sv += fx; st += us; si += iters; mn = min(mn, fx); mx = max(mx, fx);
            }
            if (check(in, bestTour) != bestVal) { cerr << "checker mismatch\n"; return 1; }
            best << name << ";" << method << ";" << bestVal;
            for (int v : bestTour) best << ";" << v;
            best << "\n";
            printf("%s %-38s avg=%9.1f min=%7lld max=%7lld  time=%7.1f ms  iters=%.1f\n", name.c_str(),
                   method.c_str(), sv / runs, mn, mx, st / runs / 1000, si / runs);
            fflush(stdout);
        }
    }
    return 0;
}
