# Evolutionary Computation - Assignment 7 (large neighborhood search)

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

Local search (LS): steepest local search with the list of improving moves from assignment 5
(inter-route node exchange + intra-route two-edges exchange).

### Large neighborhood search (LNS)

```
x ← random solution
x ← LS(x)
repeat until time limit (average time of MSLS on this instance):
    y ← destroy(x)
    y ← repair(y)
    y ← LS(y)                                 (only in the version with LS)
    if f(y) < f(x): x ← y
return x
```

### Destroy

All versions remove m = 30% of the nodes of x (m = 30 for K = 100). The remaining nodes keep their
order in the cycle.

```
destroy_random_nodes(x):
    repeat m times:
        remove a random node from x
    return x

destroy_single_subpath(x):
    s ← random position
    remove x[s], x[s+1], ..., x[s+m-1]        (positions mod K)
    return x

destroy_3_5_subpaths(x):
    p ← random from {3, 4, 5}
    r1..rp ← random split of m into p parts ≥ 1          (lengths of removed subpaths)
    k1..kp ← random split of K - m into p parts ≥ 1      (lengths of kept parts between them)
    s ← random position
    from s along the cycle: remove r1 nodes, keep k1, remove r2, keep k2, ..., remove rp, keep kp
    return x
```

The kept parts are at least 1 node long, so the removed subpaths never join into a longer one.

### Repair

Weighted 2-regret greedy cycle (the best constructive heuristic from assignment 2, weights 1:2),
started from the partial cycle after destroy:

```
repair(y):
    while |y| < K:
        for each node u not in y:
            best(u), second(u) ← two smallest  d(a,u) + d(u,b) - d(a,b)  over edges (a,b) of y
            score(u) ← (second(u) - best(u)) - 2 * (best(u) + c(u))
        u ← node with the largest score (ties: smaller best(u) + c(u))
        insert u into its best edge
    return y
```

For every free node the two best edges are remembered, after an insertion into edge (a, b) only the two
new edges are checked (full recomputation only for nodes whose remembered edge was (a, b)).

Each version run 20 times on each instance. Time limit = average running time of MSLS from assignment 6
(TSPA 712.7 ms, TSPB 691.6 ms). C++ (g++ -O2).

## Results

Cells: average (min – max). Bold = best average and best minimum in the column, underlined = methods
from the current assignment.

### Objective function

Best constructive heuristic and best local search versions from assignments 2, 3 and 5 (200 runs),
MSLS and ILS from assignment 6 (20 runs), LNS (20 runs).

{{OBJ_TABLE}}

### Running time [ms] per run

{{TIME_TABLE}}

### Number of iterations of the main loop

{{ITERS_TABLE}}

## Best solutions

Color of a node = its cost (bigger, outlined dots are the selected nodes).

{{GRIDS}}

## Checker

Yes, the best solutions were checked with the solution checker (in C++ and independently in Python from
the coordinates), values are the same.

## Code

lns.cpp (C++, standard library only). Link: https://github.com/tomasz-pres/evolutionary-computing.git

## Conclusions

- all LNS versions are much better than MSLS in the same time (TSPA 71693-72148 vs 73614, TSPB
  48685-49637 vs 50140)
- LNS is a bit worse than ILS: the best versions TSPA 71693 vs 71607, TSPB 48685 vs 48255
- 3-5 subpaths is the best destroy operator on both instances, the worst is random nodes on TSPA
  and single subpath on TSPB
- LS after repair helps only a little or not at all (TSPA 3-5 subpaths 71693 vs 71705, TSPB 3-5
  subpaths 48860 vs 48685). The repair heuristic already gives a good solution, so LS fixes little,
  and without LS there are 1.3-1.5x more iterations
- some runs get stuck (TSPB max up to 51193), the repair is deterministic, so after the same
  destroy it builds the same solution again
