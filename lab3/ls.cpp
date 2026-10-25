// Assignment 3: local search (steepest / greedy, nodes / edges intra-route moves,
// random / greedy starting solutions) + all methods from previous assignments.
//
// Build: g++ -O2 -std=c++17 -o ls ls.cpp
// Usage: ./ls [--data-dir DIR] [--out DIR] [--runs N] [--seed S]
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
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
    // distance matrix is computed once; optimisation methods use only D and cost
    in.D.assign(in.n, vector<int>(in.n, 0));
    for (int i = 0; i < in.n; i++)
        for (int j = 0; j < in.n; j++)
            in.D[i][j] = (int)floor(sqrt((in.x[i] - in.x[j]) * (in.x[i] - in.x[j]) +
                                         (in.y[i] - in.y[j]) * (in.y[i] - in.y[j])) + 0.5);
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

// ---------------------------------------------------------------- random
typedef mt19937_64 Rng;
// uniform in [0, m): multiply-shift instead of a (slow) 64-bit modulo; the bias is below 2^-32
inline int rnd(Rng &r, int m) { return (int)(((r() >> 32) * (unsigned long long)m) >> 32); }

Tour random_solution(const Instance &in, Rng &r) {
    Tour p(in.n);
    for (int i = 0; i < in.n; i++) p[i] = i;
    for (int i = 0; i < in.k; i++) swap(p[i], p[i + rnd(r, in.n - i)]);  // partial Fisher-Yates
    p.resize(in.k);
    return p;
}

// ---------------------------------------------------------------- construction heuristics
Tour nn_end(const Instance &in, int start) {
    Tour path = {start};
    vector<char> used(in.n, 0);
    used[start] = 1;
    while ((int)path.size() < in.k) {
        int best = -1, bd = INT_MAX;
        for (int v = 0; v < in.n; v++)
            if (!used[v] && in.D[path.back()][v] + in.cost[v] < bd) { bd = in.D[path.back()][v] + in.cost[v]; best = v; }
        path.push_back(best);
        used[best] = 1;
    }
    return path;
}

// Insertion cost of v into a path (rows 0..len: front, inner, back) or a cycle (rows = edges).
inline int ins_cost(const Instance &in, const Tour &s, int row, int v, bool cycle) {
    int len = s.size();
    if (cycle) {
        int a = s[row], b = s[(row + 1) % len];
        return in.D[a][v] + in.D[v][b] - in.D[a][b] + in.cost[v];
    }
    if (row == 0) return in.D[v][s[0]] + in.cost[v];
    if (row == len) return in.D[s[len - 1]][v] + in.cost[v];
    int a = s[row - 1], b = s[row];
    return in.D[a][v] + in.D[v][b] - in.D[a][b] + in.cost[v];
}

void init_seq(const Instance &in, int start, bool cycle, Tour &seq, vector<char> &used) {
    used.assign(in.n, 0);
    seq = {start};
    used[start] = 1;
    if (cycle) {
        int best = -1, bd = INT_MAX;
        for (int v = 0; v < in.n; v++)
            if (!used[v] && in.D[start][v] + in.cost[v] < bd) { bd = in.D[start][v] + in.cost[v]; best = v; }
        seq.push_back(best);
        used[best] = 1;
    }
}

// nearest-neighbour with insertion at any position / greedy cycle: cheapest insertion overall
Tour cheapest_insertion(const Instance &in, int start, bool cycle) {
    Tour seq; vector<char> used;
    init_seq(in, start, cycle, seq, used);
    while ((int)seq.size() < in.k) {
        int rows = cycle ? seq.size() : seq.size() + 1;
        int br = -1, bv = -1, bd = INT_MAX;
        for (int r = 0; r < rows; r++)
            for (int v = 0; v < in.n; v++) {
                if (used[v]) continue;
                int d = ins_cost(in, seq, r, v, cycle);
                if (d < bd) { bd = d; br = r; bv = v; }
            }
        seq.insert(seq.begin() + (cycle ? br + 1 : br), bv);
        used[bv] = 1;
    }
    return seq;
}

// choose node with the largest  wr * (2nd best insertion - best insertion) - wb * best insertion
Tour regret_construct(const Instance &in, int start, bool cycle, int wr, int wb) {
    Tour seq; vector<char> used;
    init_seq(in, start, cycle, seq, used);
    while ((int)seq.size() < in.k) {
        int rows = cycle ? seq.size() : seq.size() + 1;
        long long bscore = LLONG_MIN; int bv = -1, bpos = -1, bb1 = 0;
        for (int v = 0; v < in.n; v++) {
            if (used[v]) continue;
            int b1 = INT_MAX, b2 = INT_MAX, pos = -1;
            for (int r = 0; r < rows; r++) {
                int d = ins_cost(in, seq, r, v, cycle);
                if (d < b1) { b2 = b1; b1 = d; pos = r; }
                else if (d < b2) b2 = d;
            }
            long long score = (long long)wr * (b2 - b1) - (long long)wb * b1;
            if (score > bscore || (score == bscore && b1 < bb1)) { bscore = score; bv = v; bpos = pos; bb1 = b1; }
        }
        seq.insert(seq.begin() + (cycle ? bpos + 1 : bpos), bv);
        used[bv] = 1;
    }
    return seq;
}

// ---------------------------------------------------------------- local search
//
// Neighbourhood = inter-route moves (swap selected node at tour position i with an unselected
// node) + intra-route moves (either swap of two nodes, or exchange of two edges = reversal of
// the segment between them). All moves are indexed in one combined space [0, M):
//   [0, K*(n-K))                  inter moves  (position i, unselected slot j)
//   [K*(n-K), M)                  intra moves  (index into a precomputed list of valid pairs)
// Only the delta of the objective is computed for a move, never the whole objective.

struct LocalSearch {
    const Instance &in;
    bool steepest, edges;
    int K, U;                     // selected / unselected counts
    Tour t;                       // current tour
    Tour unsel;                   // unselected nodes
    vector<int> pa, pb;           // valid intra pairs
    long long M, nInter;
    long long evals = 0;

    LocalSearch(const Instance &in_, bool steepest_, bool edges_) : in(in_), steepest(steepest_), edges(edges_) {
        K = in.k; U = in.n - K;
        for (int i = 0; i < K; i++)
            for (int j = i + 1; j < K; j++) {
                if (edges && (j < i + 2 || (i == 0 && j == K - 1))) continue;  // edges sharing a node: no-op
                pa.push_back(i); pb.push_back(j);
            }
        nInter = (long long)K * U;
        M = nInter + pa.size();
    }

    inline int at(int i) const { return t[(i % K + K) % K]; }

    // delta of replacing t[i] by unselected node v
    inline int delta_inter(int i, int v) const {
        int a = t[i], p = at(i - 1), nx = at(i + 1);
        return in.D[p][v] + in.D[v][nx] - in.D[p][a] - in.D[a][nx] + in.cost[v] - in.cost[a];
    }
    // delta of swapping the nodes at positions i < j
    inline int delta_swap(int i, int j) const {
        if (i == 0 && j == K - 1) { i = K - 1; j = 0; }       // cyclically adjacent
        int a = t[i], b = t[j];
        int pi = at(i - 1), nj = at(j + 1);
        if ((i + 1) % K == j)                                   // adjacent: ... pi a b nj ...
            return in.D[pi][b] + in.D[a][nj] - in.D[pi][a] - in.D[b][nj];
        int ni = at(i + 1), pj = at(j - 1);
        return in.D[pi][b] + in.D[b][ni] + in.D[pj][a] + in.D[a][nj]
             - in.D[pi][a] - in.D[a][ni] - in.D[pj][b] - in.D[b][nj];
    }
    // delta of exchanging edges (i,i+1) and (j,j+1): new edges (i,j), (i+1,j+1)
    inline int delta_edges(int i, int j) const {
        int a = t[i], a1 = t[i + 1], b = t[j], b1 = at(j + 1);
        return in.D[a][b] + in.D[a1][b1] - in.D[a][a1] - in.D[b][b1];
    }

    inline int delta(long long m) {
        evals++;
        if (m < nInter) return delta_inter((int)(m / U), unsel[m % U]);
        int idx = (int)(m - nInter);
        return edges ? delta_edges(pa[idx], pb[idx]) : delta_swap(pa[idx], pb[idx]);
    }
    inline void apply(long long m) {
        if (m < nInter) {
            int i = (int)(m / U), j = (int)(m % U);
            swap(t[i], unsel[j]);
            return;
        }
        int idx = (int)(m - nInter);
        if (edges) reverse(t.begin() + pa[idx] + 1, t.begin() + pb[idx] + 1);
        else swap(t[pa[idx]], t[pb[idx]]);
    }

    // returns final tour; evals counts the evaluated moves
    Tour run(const Tour &start, Rng &rng) {
        t = start; evals = 0;
        vector<char> sel(in.n, 0);
        for (int v : t) sel[v] = 1;
        unsel.clear();
        for (int v = 0; v < in.n; v++) if (!sel[v]) unsel.push_back(v);

        if (steepest) {
            while (true) {
                long long bm = -1; int bd = 0;
                for (long long m = 0; m < M; m++) {
                    int d = delta(m);
                    if (d < bd) { bd = d; bm = m; }
                }
                if (bm < 0) break;
                apply(bm);
            }
        } else {
            // greedy: every iteration browses the whole (inter + intra) move index space in a
            // fresh uniformly random order. The permutation is generated lazily (Fisher-Yates,
            // one step per examined move), so the cost is proportional to the number of
            // evaluated moves. Deltas are computed one move at a time and the first improving
            // move is applied.
            vector<int> perm(M);
            for (long long m = 0; m < M; m++) perm[m] = (int)m;
            bool improved = true;
            while (improved) {
                improved = false;
                for (long long q = 0; q < M; q++) {
                    swap(perm[q], perm[q + rnd(rng, (int)(M - q))]);
                    if (delta(perm[q]) < 0) { apply(perm[q]); improved = true; break; }
                }
            }
        }
        return t;
    }
};

// ---------------------------------------------------------------- experiment
struct Method {
    string name;
    bool isNew;
    // produces a solution for run r; evals filled for local-search methods
    function<Tour(int r, long long &evals)> fn;
};

struct RunResult { long long value; double us; long long evals; };
double g_ls_us = 0;  // time of the local search part of the last run (without the starting solution)

int main(int argc, char **argv) {
    string dataDir = "../lab1", outDir = "results";
    int runs = 200; unsigned long long seed = 42;
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        if (a == "--data-dir" && i + 1 < argc) dataDir = argv[++i];
        else if (a == "--out" && i + 1 < argc) outDir = argv[++i];
        else if (a == "--runs" && i + 1 < argc) runs = atoi(argv[++i]);
        else if (a == "--seed" && i + 1 < argc) seed = strtoull(argv[++i], 0, 10);
        else { cerr << "usage: ls [--data-dir DIR] [--out DIR] [--runs N] [--seed S]\n"; return 1; }
    }
    string cmd = "mkdir -p '" + outDir + "'";
    if (system(cmd.c_str()) != 0) return 1;
    ofstream raw(outDir + "/raw.csv"), best(outDir + "/best.txt");
    raw << "instance,method,run,value,time_us,evals,ls_time_us\n";

    for (string name : {"TSPA", "TSPB"}) {
        Instance in = load_instance(dataDir, name);
        const int n = in.n;

        // the same random starting solutions are used by every method that needs them
        auto random_start = [&](int r) { Rng g(seed * 1000003ULL + r); return random_solution(in, g); };
        auto greedy_start = [&](int s) { return regret_construct(in, s % n, true, 1, 2); };

        vector<Method> methods;
        methods.push_back({"Random", false, [&](int r, long long &) { return random_start(r); }});
        methods.push_back({"NN (end only)", false, [&](int r, long long &) { return nn_end(in, r % n); }});
        methods.push_back({"NN (any position)", false, [&](int r, long long &) { return cheapest_insertion(in, r % n, false); }});
        methods.push_back({"Greedy cycle", false, [&](int r, long long &) { return cheapest_insertion(in, r % n, true); }});
        methods.push_back({"2-regret (NN any)", false, [&](int r, long long &) { return regret_construct(in, r % n, false, 1, 0); }});
        methods.push_back({"2-regret (greedy cycle)", false, [&](int r, long long &) { return regret_construct(in, r % n, true, 1, 0); }});
        methods.push_back({"Weighted (NN any)", false, [&](int r, long long &) { return regret_construct(in, r % n, false, 1, 1); }});
        methods.push_back({"Weighted (greedy cycle)", false, [&](int r, long long &) { return regret_construct(in, r % n, true, 1, 1); }});

        methods.push_back({"Weighted (greedy cycle, 1:2)", false, [&](int r, long long &) { return regret_construct(in, r % n, true, 1, 2); }});

        // local search objects are reused between runs (neighbourhood tables built once, untimed)
        for (int steepest = 1; steepest >= 0; steepest--)
            for (int edges = 0; edges <= 1; edges++)
                for (int greedyStart = 0; greedyStart <= 1; greedyStart++) {
                    auto ls = make_shared<LocalSearch>(in, steepest, edges);
                    string nm = string(steepest ? "Steepest" : "Greedy") + (edges ? ", edges" : ", nodes") +
                                (greedyStart ? ", greedy start" : ", random start");
                    int sd = steepest * 4 + edges * 2 + greedyStart;
                    methods.push_back({nm, true, [&in, ls, greedyStart, random_start, greedy_start, seed, sd](int r, long long &ev) {
                        Tour s = greedyStart ? greedy_start(r) : random_start(r);
                        Rng g(seed * 7919ULL + sd * 100003ULL + r);
                        auto c0 = chrono::steady_clock::now();
                        Tour res = ls->run(s, g);
                        g_ls_us = chrono::duration<double, micro>(chrono::steady_clock::now() - c0).count();
                        ev = ls->evals;
                        return res;
                    }});
                }

        for (auto &m : methods) {
            vector<RunResult> res;
            Tour bestTour; long long bestVal = LLONG_MAX;
            for (int r = 0; r < runs; r++) {
                long long ev = 0;
                g_ls_us = 0;
                auto t0 = chrono::steady_clock::now();
                Tour t = m.fn(r, ev);
                auto t1 = chrono::steady_clock::now();
                double us = chrono::duration<double, micro>(t1 - t0).count();
                long long v = objective(in, t);
                res.push_back({v, us, ev});
                raw << name << ",\"" << m.name << "\"," << r << "," << v << "," << us << "," << ev << "," << g_ls_us << "\n";
                if (v < bestVal) { bestVal = v; bestTour = t; }
            }
            if (check(in, bestTour) != bestVal) { cerr << "checker mismatch\n"; return 1; }
            best << name << ";" << m.name << ";" << bestVal;
            for (int v : bestTour) best << ";" << v;
            best << "\n";
            double sv = 0, st = 0;
            long long mn = LLONG_MAX, mx = 0;
            for (auto &x : res) { sv += x.value; st += x.us; mn = min(mn, x.value); mx = max(mx, x.value); }
            printf("%s %-40s avg=%9.1f min=%7lld max=%7lld  time=%9.1f us\n", name.c_str(), m.name.c_str(),
                   sv / runs, mn, mx, st / runs);
            fflush(stdout);
        }
    }
    return 0;
}
