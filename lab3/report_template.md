# Evolutionary Computation - Assignment 3 (local search)

Authors:
1. Name: Tomasz Preś  Index number: 160305
2. Name: Marek Olejniczak  Index number: 160303

## Problem

We are given three columns of integers with a row for each node. The first two columns contain x
and y coordinates of the node positions in a plane. The third column contains node costs. The goal is
to select exactly 50% of the nodes (if the number of nodes is odd we round the number of nodes to
be selected up) and form a Hamiltonian cycle (closed path) through this set of nodes such that the
sum of the total length of the path plus the total cost of the selected nodes is minimized.

The distances between nodes are calculated as Euclidean distances rounded mathematically to
integer values. The distance matrix is calculated just after reading an instance and then only
the distance matrix (no nodes coordinates) is accessed by the optimization methods.

## Algorithms (pseudocode)

All methods from assignments 1 and 2 are included in the tables (re-implemented in C++, so that all
running times were measured in the same way; the objective values of the constructive methods are
identical to the previous reports, only *Random* differs because another random generator was used).

### Neighborhood

The solution is a cycle `t[0..K-1]` (K = 100) of selected nodes; the other nodes are unselected.
Three kinds of moves, only the change of the objective (delta) is calculated, never the whole objective.

```
d = distance, c = node cost, all indices modulo K.
For a node x at position i:  pred(x) = t[i-1],  succ(x) = t[i+1].

inter-route (K*(n-K) = 10000 moves):
    replace selected a = t[i] by an unselected v
    delta = d(pred(a), v) + d(v, succ(a))
          - d(pred(a), a) - d(a, succ(a))
          + c(v) - c(a)

intra-route, two nodes exchange (4950 moves):
    swap a = t[i] and b = t[j]
    delta = d(pred(a), b) + d(b, succ(a))
          + d(pred(b), a) + d(a, succ(b))
          - d(pred(a), a) - d(a, succ(a))
          - d(pred(b), b) - d(b, succ(b))
    (if a and b are neighbours, the edge between them
     is unchanged and is not counted)

intra-route, two edges exchange (4850 moves):
    remove edges (t[i], t[i+1]) and (t[j], t[j+1]),
    reconnect by reversing the segment t[i+1..j]
    delta = d(t[i], t[j]) + d(t[i+1], t[j+1])
          - d(t[i], t[i+1]) - d(t[j], t[j+1])
```

(Pairs of edges that share a node give no change of the cycle, so they are excluded from the edge exchange.)
The inter-route moves and the intra-route moves together are one neighborhood of ~14900 moves.
Cost of a node selected / not selected changes only in the inter-route move.

### Steepest local search

```
x = starting solution
loop:
    best = none
    for each move m in the whole neighborhood (inter + intra):
        if delta(m) < best.delta (and < 0): best = m
    if best is none: stop (local optimum)
    apply best
```

### Greedy local search

```
x = starting solution
loop:
    browse all moves (inter + intra together) in a random order:
        evaluate delta(m) of the next move
        if delta(m) < 0: apply m, go to the next iteration of the loop
    if no move improved: stop (local optimum)
```

**Randomization.** The moves of both kinds are indexed by one common range `0 .. M-1`
(first the 10000 inter-route moves, then the intra-route pairs). In every iteration the
index range is browsed in a new, uniformly random permutation, so the two kinds of moves are
mixed (it never starts from a given kind). The permutation is not generated in advance: it is
produced lazily with the Fisher-Yates shuffle, one step for every move that is going to be
evaluated (`swap(perm[q], perm[q + rand(M - q)])`, then evaluate `perm[q]`). Thanks to this,
only the moves up to the first improving one are evaluated, and also the cost of the randomization
itself is proportional to the number of evaluated moves (our first version shuffled the whole
array in every iteration, which cost as much as evaluating all deltas, so greedy was not faster).

### Starting solutions

- random: 200 random solutions (random 100 nodes in random order); the same 200 solutions are used
  by all four methods with random starts, so the comparison is paired,
- greedy: the best method from assignments 1-2, i.e. the weighted sum regret heuristic (greedy cycle,
  weights 1:2, i.e. regret - 2 * best insertion cost, which had the best average on both instances), started from each of the 200 nodes.

Every method was run 200 times. The time of one run includes creating the starting solution
(for greedy starts this is the construction heuristic, about 1.5 ms of the total). The additional
table shows the number of evaluated moves (deltas) in one run, which does not depend on the
implementation. C++ (g++ -O2), one thread, wall-clock time of every run measured separately.

## Results

Cells: average (min – max). Bold = best average and best minimum in the column, underlined = methods
from the current assignment.

### Objective function

{{OBJ_TABLE}}

### Running time [ms] per run

{{TIME_TABLE}}

### Running time [ms] per run, local search only (without creating the starting solution)

{{LSTIME_TABLE}}

### Number of evaluated moves per run (local search only)

{{EVALS_TABLE}}

## Best solutions

Color of a node = its cost (bigger, outlined dots are the selected nodes).

{{GRIDS}}

## Checker

Yes, all the best solutions were checked with the solution checker (in C++ and independently in
Python from the coordinates), values are the same.

## Code

ls.cpp (C++, standard library only). Link: https://github.com/tomasz-pres/evolutionary-computing.git

## Conclusions

- **Edges vs nodes.** The two-edges exchange is much better from random starts (avg TSPA ~76.2k vs ~88-90k,
  TSPB ~52.6-52.9k vs ~63.8-66.2k): reversing a segment removes the crossing edges that random tours are full
  of. From greedy starts the difference is negligible (TSPA 72811 vs 72862).
- **Steepest vs greedy: quality.** Practically the same (differences below 1% from random starts, ~0.1% from
  greedy starts), so taking the first improving move instead of the best one costs almost nothing.
- **Steepest vs greedy: time.** From random starts greedy is about 9x faster (TSPA nodes 2.3 vs 21.1 ms, edges
  1.4 vs 12.9 ms) and evaluates 17-20x fewer moves; the time gain is smaller because one greedy evaluation is
  ~2x more expensive (random number, random memory access). From greedy starts the search is short
  (0.4-1.2 ms of ~2-2.7 ms) and dominated by the construction (~1.5 ms); both versions also have to scan the
  whole neighborhood once at the end to confirm the local optimum, so greedy saves only ~2x evaluations, which
  the costlier evaluation cancels. Greedy is therefore faster in 3 of 4 such cases, but not in TSPA edges
  (1.98 vs 1.91 ms).
- **Random vs greedy starts.** Greedy starts are much better on TSPA and faster to improve
  (~2 ms vs 13-21 ms). The exception is TSPB with edges: random starts (avg 52.6k steepest, 52.9k greedy)
  beat greedy starts (54.8-55.1k) and all constructive methods (best 55.7k), and give the best TSPB solution
  (49873) - the greedy start sticks near its own local optimum, while random starts are more diverse.
- **Gain over the heuristic.** Local search from the greedy start improves it only a little (TSPA 73191 -> 72820,
  TSPB 55668 -> 54788) but cheaply (+0.4-1.2 ms).
- **Overall.** Best average on TSPA: steepest with edges from greedy starts (72811, greedy 72820); on TSPB:
  steepest with edges from random starts (52576). Greedy with edges from random starts is only 0.7% worse
  on TSPB and 9x faster, so it has the best quality/time trade-off.
