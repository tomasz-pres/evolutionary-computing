# Evolutionary Computation - Assignment 4 (candidate moves in local search)

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

We use the neighborhood that was the best in assignment 3: inter-route exchange (selected node <->
not selected node) + intra-route two-edges exchange, steepest version, random starting solutions.

### Baseline: steepest local search (from assignment 3)

```
x = random solution
loop:
    evaluate ALL moves: K*(n-K) = 10000 node exchanges + 4850 two-edges exchanges
    if the best delta < 0: apply it, else stop
```

Same code as in assignment 3 (only the distance matrix is stored as one flat array), it gives exactly
the same results as there (same random starting solutions, same values and numbers of evaluated moves
in all 400 runs).

### Candidate edges

```
for each node a:
    near[a] = the 10 nodes b != a with the smallest  d(a, b) + cost(b)
candidate edge = (a, b) with b in near[a]
```

The lists are computed once per instance, before the local search (not counted in the time).

### Candidate moves

A move is a candidate move if it adds at least one candidate edge to the solution. We do not go through
all moves and filter them, we generate only the candidate moves:

```
for each selected node a (position i):        pred/succ = neighbours in the cycle
    for each b in near[a]:
        if b is selected:
            # two-edges exchanges that add edge (a, b)
            remove (a, succ a), (b, succ b), add (a, b), (succ a, succ b)
                delta = d(a,b) + d(succ a, succ b) - d(a, succ a) - d(b, succ b)
            remove (pred a, a), (pred b, b), add (pred a, pred b), (a, b)
                delta = d(pred a, pred b) + d(a,b) - d(pred a, a) - d(pred b, b)
            (skipped if the two removed edges share a node)
        else:
            # node exchanges that add edge (a, b) - b replaces a NEIGHBOUR of a, not a itself
            replace succ a by b: new edges (a, b), (b, succ succ a)
                delta = d(a,b) + d(b, succ succ a) - d(a, succ a) - d(succ a, succ succ a) + c(b) - c(succ a)
            replace pred a by b: new edges (pred pred a, b), (b, a)
                delta = d(pred pred a, b) + d(b,a) - d(pred pred a, pred a) - d(pred a, a) + c(b) - c(pred a)
```

So there are 2 moves for every candidate edge, ~K * 10 * 2 = 2000 moves per iteration instead of 14850.
Positions of the nodes are kept in an array `pos[node]` (-1 = not selected), so checking if b is
selected and finding its neighbours is O(1).

### Steepest local search with candidate moves

```
x = random solution
loop:
    evaluate all candidate moves (above)
    if the best delta < 0: apply it (update pos[]), else stop
```

The size of the lists (10) was also tested experimentally: 5, 10, 15, 20.

**Correctness checks.** (1) with lists of all 199 nodes the candidate moves are the whole neighborhood:
all 400 results were local optima of the full neighborhood (results differ from the baseline only in
23/200 runs because of ties - equal deltas evaluated in a different order); (2) in a test build (50 runs of every method) every
applied move was checked: the real change of the objective = the computed delta.

Every method was run 200 times from the same 200 random solutions (the same ones as in assignment 3).
C++ (g++ -O2), one thread, wall-clock time of every run measured separately.

## Results

Cells: average (min – max). Bold = best average and best minimum in the column, underlined = methods
from the current assignment.

### Objective function

Methods from assignments 1-3: values from assignment 3 (they do not depend on the computer).

{{OBJ_TABLE}}

### Running time [ms] per run

Only the methods run in this assignment (on the same computer).

{{TIME_TABLE}}

### Running time [ms] per run, local search only (without creating the starting solution)

{{LSTIME_TABLE}}

### Number of evaluated moves per run

{{EVALS_TABLE}}

## Best solutions

Color of a node = its cost (bigger, outlined dots are the selected nodes).

{{GRIDS}}

## Checker

Yes, the best solutions were checked with the solution checker (in C++ and independently in Python from
the coordinates), values are the same.

## Code

cand.cpp (C++, standard library only). Link: https://github.com/tomasz-pres/evolutionary-computing.git

## Conclusions

- candidate moves with 10 neighbours make the local search about 5x faster (TSPA 3.21 -> 0.63 ms,
  TSPB 3.20 -> 0.64 ms), they evaluate about 8x fewer moves (~250k instead of ~2M). One candidate move
  is a bit more expensive to evaluate (position lookups), so the time gain is smaller than the gain in moves
- on TSPB the quality is the same as without candidates (avg 52565 vs 52576)
- on TSPA the quality is worse (avg 80057 vs 76237, +5%). TSPA has bigger node costs (0-2000 vs 0-1000 in
  TSPB), so "nearest" with cost is less geometric: only 36% of the 10 candidates are among the 10
  geometrically closest nodes (50% in TSPB). Two-edges exchange does not depend on node costs, so
  some good short edges are never tried
- bigger lists fix it but are slower: TSPA with 15 -> 77511 (3.4x faster), with 20 -> 76607
  (+0.5%, 2.3x faster). With 5 the search is ~10x faster but clearly worse on both instances
- so 10 is a good choice for TSPB, for TSPA 15-20 is a better trade-off
