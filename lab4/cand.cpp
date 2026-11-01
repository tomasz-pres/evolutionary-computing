// Assignment 4: candidate moves in steepest local search.
// Neighbourhood = the best one from assignment 3: inter-route node exchange (selected <-> unselected)
// + intra-route two-edges exchange. Baseline = steepest over the whole neighbourhood (as in lab3),
// new method = steepest over candidate moves only (moves introducing at least one candidate edge).
//
// Build: g++ -O2 -std=c++17 -o cand cand.cpp
// Usage: ./cand [--data-dir DIR] [--out DIR] [--runs N] [--seed S]
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
typedef vector<int> Tour;

// ---------------------------------------------------------------- instance (same as lab3)
struct Instance {
    string name;
    int n = 0, k = 0;
    vector<double> x, y;
    vector<int> cost;
    vector<vector<int>> D;
    vector<int> F;   // the same distances as one flat array: F[a * n + b] (faster access)
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

// same generator and seeds as lab3, so the starting solutions are the same as there
Tour random_solution(const Instance &in, Rng &r) {
    Tour p(in.n);
    for (int i = 0; i < in.n; i++) p[i] = i;
    for (int i = 0; i < in.k; i++) swap(p[i], p[i + rnd(r, in.n - i)]);
    p.resize(in.k);
    return p;
}

// ---------------------------------------------------------------- baseline: full steepest (lab3)
struct FullSteepest {
    const Instance &in;
    int K, U;
    Tour t, unsel;
    vector<int> pa, pb;   // valid two-edges exchange pairs (i < j, edges not sharing a node)
    long long evals = 0;

    FullSteepest(const Instance &in_) : in(in_) {
        K = in.k; U = in.n - K;
        for (int i = 0; i < K; i++)
            for (int j = i + 2; j < K; j++)
                if (!(i == 0 && j == K - 1)) { pa.push_back(i); pb.push_back(j); }
    }
    inline int at(int i) const { return t[(i % K + K) % K]; }

    Tour run(const Tour &start) {
        t = start; evals = 0;
        vector<char> sel(in.n, 0);
        for (int v : t) sel[v] = 1;
        unsel.clear();
        for (int v = 0; v < in.n; v++) if (!sel[v]) unsel.push_back(v);
        const int n = in.n; const int *F = in.F.data(); const int *C = in.cost.data();
        while (true) {
            int bd = 0, kind = -1, bi = 0, bj = 0;
            // inter-route: replace t[i] by unselected unsel[j]
            for (int i = 0; i < K; i++) {
                int a = t[i], p = at(i - 1), nx = at(i + 1);
                const int *Fp = F + p * n, *Fn = F + nx * n;
                int base = -Fp[a] - Fn[a] - C[a];
                for (int j = 0; j < U; j++) {
                    int v = unsel[j];
                    int d = base + Fp[v] + Fn[v] + C[v];
                    if (d < bd) { bd = d; kind = 0; bi = i; bj = j; }
                }
                evals += U;
            }
            // intra-route: exchange edges (i,i+1), (j,j+1)
            for (size_t q = 0; q < pa.size(); q++) {
                int i = pa[q], j = pb[q];
                int a = t[i], a1 = t[i + 1], b = t[j], b1 = at(j + 1);
                int d = F[a * n + b] + F[a1 * n + b1] - F[a * n + a1] - F[b * n + b1];
                if (d < bd) { bd = d; kind = 1; bi = i; bj = j; }
            }
            evals += pa.size();
            if (kind < 0) break;
            if (kind == 0) swap(t[bi], unsel[bj]);
            else reverse(t.begin() + bi + 1, t.begin() + bj + 1);
        }
        return t;
    }
};

// ---------------------------------------------------------------- candidate moves steepest
//
// Candidate edges: (a, b) for b in near[a], where near[a] = the `nc` nodes with the smallest
// D[a][b] + cost[b]. A move is a candidate move if it introduces a candidate edge into the solution.
// For every selected node a and every b in near[a]:
//   b selected   -> two two-edges exchanges introducing edge (a, b):
//                   remove (a, succ a), (b, succ b), add (a, b), (succ a, succ b)
//                   remove (pred a, a), (pred b, b), add (pred a, pred b), (a, b)
//   b unselected -> two node exchanges introducing edge (a, b):
//                   replace succ(a) by b  -> new edges (a, b), (b, succ succ a)
//                   replace pred(a) by b  -> new edges (pred pred a, b), (b, a)
// Only these moves are evaluated; the full neighbourhood is never browsed.
struct CandidateSteepest {
    const Instance &in;
    int K;
    vector<vector<int>> near;
    Tour t;
    vector<int> pos;      // position in t, -1 if not selected
    long long evals = 0;

    CandidateSteepest(const Instance &in_, int nc) : in(in_) {
        K = in.k;
        near.assign(in.n, {});
        for (int a = 0; a < in.n; a++) {
            vector<int> o;
            for (int b = 0; b < in.n; b++) if (b != a) o.push_back(b);
            stable_sort(o.begin(), o.end(), [&](int u, int v) {
                return in.D[a][u] + in.cost[u] < in.D[a][v] + in.cost[v];
            });
            o.resize(min<int>(nc, o.size()));
            near[a] = o;
        }
    }
    Tour run(const Tour &start) {
        t = start; evals = 0;
        pos.assign(in.n, -1);
        for (int i = 0; i < K; i++) pos[t[i]] = i;
        const int n = in.n; const int *F = in.F.data(); const int *C = in.cost.data();
        while (true) {
            // kind 0: replace node at position x by unselected node y; kind 1: two-edges exchange of
            // edges starting at positions x and y
            int bd = 0, kind = -1, bx = 0, by = 0;
            for (int i = 0; i < K; i++) {
                int ip = i ? i - 1 : K - 1, in_ = i + 1 < K ? i + 1 : 0;
                int a = t[i], pa = t[ip], na = t[in_];
                const int *Fa = F + a * n;
                int d_a_na = Fa[na], d_pa_a = Fa[pa];
                for (int b : near[a]) {
                    int j = pos[b];
                    const int *Fb = F + b * n;
                    if (j >= 0) {
                        int jp = j ? j - 1 : K - 1, jn = j + 1 < K ? j + 1 : 0;
                        int pb = t[jp], nb = t[jn];
                        // remove (a, na), (b, nb), add (a, b), (na, nb)
                        if (b != na && nb != a) {
                            int d = Fa[b] + F[na * n + nb] - d_a_na - Fb[nb];
                            if (d < bd) { bd = d; kind = 1; bx = i; by = j; }
                        }
                        // remove (pa, a), (pb, b), add (pa, pb), (a, b)
                        if (pb != a && pa != b) {
                            int d = F[pa * n + pb] + Fa[b] - d_pa_a - Fb[pb];
                            if (d < bd) { bd = d; kind = 1; bx = ip; by = jp; }
                        }
                    } else {
                        // replace na by b: new edges (a, b), (b, nna)
                        int nna = t[in_ + 1 < K ? in_ + 1 : 0];
                        int d = Fa[b] + Fb[nna] - d_a_na - F[na * n + nna] + C[b] - C[na];
                        if (d < bd) { bd = d; kind = 0; bx = in_; by = b; }
                        // replace pa by b: new edges (ppa, b), (b, a)
                        int ppa = t[ip ? ip - 1 : K - 1];
                        d = F[ppa * n + b] + Fb[a] - F[ppa * n + pa] - d_pa_a + C[b] - C[pa];
                        if (d < bd) { bd = d; kind = 0; bx = ip; by = b; }
                    }
                    evals += 2;
                }
            }
            if (kind < 0) break;
#ifdef CHECK_DELTA
            long long before = objective(in, t);
#endif
            if (kind == 0) {
                int old = t[bx];
                t[bx] = by; pos[by] = bx; pos[old] = -1;
            } else {
                int i = min(bx, by), j = max(bx, by);
                reverse(t.begin() + i + 1, t.begin() + j + 1);
                for (int q = i + 1; q <= j; q++) pos[t[q]] = q;
            }
#ifdef CHECK_DELTA
            if (objective(in, t) - before != bd) { cerr << "delta mismatch\n"; exit(1); }
#endif
        }
        return t;
    }
};

// ---------------------------------------------------------------- experiment
int main(int argc, char **argv) {
    string dataDir = "../data", outDir = "results";
    int runs = 200; unsigned long long seed = 42;
    vector<int> sizes = {5, 10, 15, 20};
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        if (a == "--data-dir" && i + 1 < argc) dataDir = argv[++i];
        else if (a == "--out" && i + 1 < argc) outDir = argv[++i];
        else if (a == "--runs" && i + 1 < argc) runs = atoi(argv[++i]);
        else if (a == "--seed" && i + 1 < argc) seed = strtoull(argv[++i], 0, 10);
        else { cerr << "usage: cand [--data-dir DIR] [--out DIR] [--runs N] [--seed S]\n"; return 1; }
    }
    string cmd = "mkdir -p '" + outDir + "'";
    if (system(cmd.c_str()) != 0) return 1;
    ofstream raw(outDir + "/raw.csv"), best(outDir + "/best.txt");
    raw << "instance,method,run,value,time_us,evals,ls_time_us\n";

    for (string name : {"TSPA", "TSPB"}) {
        Instance in = load_instance(dataDir, name);
        auto random_start = [&](int r) { Rng g(seed * 1000003ULL + r); return random_solution(in, g); };

        FullSteepest full(in);
        vector<CandidateSteepest> cands;
        for (int nc : sizes) cands.emplace_back(in, nc);   // candidate lists built once, untimed

        vector<string> names = {"Steepest, edges, random start"};
        for (int nc : sizes) names.push_back("Candidates (" + to_string(nc) + "), steepest, edges, random start");

        for (size_t mi = 0; mi < names.size(); mi++) {
            Tour bestTour; long long bestVal = LLONG_MAX;
            double sv = 0, st = 0; long long mn = LLONG_MAX, mx = 0;
            for (int r = 0; r < runs; r++) {
                auto t0 = chrono::steady_clock::now();
                Tour s = random_start(r);
                auto t1 = chrono::steady_clock::now();
                Tour res; long long ev;
                if (mi == 0) { res = full.run(s); ev = full.evals; }
                else { res = cands[mi - 1].run(s); ev = cands[mi - 1].evals; }
                auto t2 = chrono::steady_clock::now();
                double us = chrono::duration<double, micro>(t2 - t0).count();
                double ls_us = chrono::duration<double, micro>(t2 - t1).count();
                long long v = objective(in, res);
                raw << name << ",\"" << names[mi] << "\"," << r << "," << v << "," << us << "," << ev << "," << ls_us << "\n";
                if (v < bestVal) { bestVal = v; bestTour = res; }
                sv += v; st += us; mn = min(mn, v); mx = max(mx, v);
            }
            if (check(in, bestTour) != bestVal) { cerr << "checker mismatch\n"; return 1; }
            best << name << ";" << names[mi] << ";" << bestVal;
            for (int v : bestTour) best << ";" << v;
            best << "\n";
            printf("%s %-52s avg=%9.1f min=%7lld max=%7lld  time=%8.3f ms\n", name.c_str(), names[mi].c_str(),
                   sv / runs, mn, mx, st / runs / 1000);
            fflush(stdout);
        }
    }
    return 0;
}
