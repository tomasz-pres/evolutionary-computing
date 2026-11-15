// MSLS and ILS with steepest local search (list of moves from assignment 5)
// Build: g++ -O2 -std=c++17 -o ils ils.cpp
// Usage: ./ils [--data-dir DIR] [--out DIR] [--runs N] [--msls-iters N] [--seed S] [--db N] [--nr N]
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
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

// ---------------------------------------------------------------- perturbation
// db double bridges: cycle cut into 4 parts A B C D -> A C B D
// nr node replacements: random node of the cycle replaced by a random node from outside
Tour perturb(const Instance &in, const Tour &x, Rng &r, int db, int nr) {
    int K = x.size();
    Tour t = x;
    for (int q = 0; q < db; q++) {
        int c[3];
        do {
            for (int &v : c) v = 1 + rnd(r, K - 1);
            sort(c, c + 3);
        } while (c[0] == c[1] || c[1] == c[2]);
        Tour y(t.begin(), t.begin() + c[0]);
        y.insert(y.end(), t.begin() + c[1], t.begin() + c[2]);
        y.insert(y.end(), t.begin() + c[0], t.begin() + c[1]);
        y.insert(y.end(), t.begin() + c[2], t.end());
        t = y;
    }
    vector<char> sel(in.n, 0);
    for (int v : t) sel[v] = 1;
    vector<int> out;
    for (int v = 0; v < in.n; v++) if (!sel[v]) out.push_back(v);
    for (int q = 0; q < nr; q++) {
        int i = rnd(r, K), j = rnd(r, out.size());
        swap(t[i], out[j]);
    }
    return t;
}

// ---------------------------------------------------------------- experiment
int main(int argc, char **argv) {
    string dataDir = "../data", outDir = "results";
    int runs = 20, msls_iters = 400, db = 1, nr = 8;
    unsigned long long seed = 42;
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        if (a == "--data-dir" && i + 1 < argc) dataDir = argv[++i];
        else if (a == "--out" && i + 1 < argc) outDir = argv[++i];
        else if (a == "--runs" && i + 1 < argc) runs = atoi(argv[++i]);
        else if (a == "--msls-iters" && i + 1 < argc) msls_iters = atoi(argv[++i]);
        else if (a == "--seed" && i + 1 < argc) seed = strtoull(argv[++i], 0, 10);
        else if (a == "--db" && i + 1 < argc) db = atoi(argv[++i]);
        else if (a == "--nr" && i + 1 < argc) nr = atoi(argv[++i]);
        else { cerr << "usage: ils [--data-dir DIR] [--out DIR] [--runs N] [--msls-iters N] [--seed S] [--db N] [--nr N]\n"; return 1; }
    }
    string cmd = "mkdir -p '" + outDir + "'";
    if (system(cmd.c_str()) != 0) return 1;
    ofstream raw(outDir + "/raw.csv"), best(outDir + "/best.txt");
    raw << "instance,method,run,value,time_us,ls_runs\n";
    typedef chrono::steady_clock Clock;

    for (string name : {"TSPA", "TSPB"}) {
        Instance in = load_instance(dataDir, name);
        LMSteepest ls(in);
        vector<string> names = {"MSLS", "ILS"};
        double msls_avg_us = 0;
        for (int mi = 0; mi < 2; mi++) {
            Tour bestTour; long long bestVal = LLONG_MAX;
            double sv = 0, st = 0, sl = 0; long long mn = LLONG_MAX, mx = 0;
            for (int r = 0; r < runs; r++) {
                Rng g(seed * 1000003ULL + r * 7919ULL + mi);
                auto t0 = Clock::now();
                Tour x; long long fx = LLONG_MAX; int lsRuns = 0;
                if (mi == 0) {   // MSLS: best of msls_iters local searches from random solutions
                    for (int it = 0; it < msls_iters; it++) {
                        Tour y = ls.run(random_solution(in, g));
                        long long fy = objective(in, y);
                        lsRuns++;
                        if (fy < fx) { fx = fy; x = y; }
                    }
                } else {         // ILS: until the average time of MSLS
                    x = ls.run(random_solution(in, g));
                    fx = objective(in, x);
                    lsRuns = 1;
                    while (chrono::duration<double, micro>(Clock::now() - t0).count() < msls_avg_us) {
                        Tour y = ls.run(perturb(in, x, g, db, nr));
                        long long fy = objective(in, y);
                        lsRuns++;
                        if (fy < fx) { fx = fy; x = y; }
                    }
                }
                double us = chrono::duration<double, micro>(Clock::now() - t0).count();
                raw << name << ",\"" << names[mi] << "\"," << r << "," << fx << "," << us << "," << lsRuns << "\n";
                if (fx < bestVal) { bestVal = fx; bestTour = x; }
                sv += fx; st += us; sl += lsRuns; mn = min(mn, fx); mx = max(mx, fx);
            }
            if (mi == 0) msls_avg_us = st / runs;
            if (check(in, bestTour) != bestVal) { cerr << "checker mismatch\n"; return 1; }
            best << name << ";" << names[mi] << ";" << bestVal;
            for (int v : bestTour) best << ";" << v;
            best << "\n";
            printf("%s %-5s avg=%9.1f min=%7lld max=%7lld  time=%8.1f ms  ls_runs=%.1f\n", name.c_str(),
                   names[mi].c_str(), sv / runs, mn, mx, st / runs / 1000, sl / runs);
            fflush(stdout);
        }
    }
    return 0;
}
