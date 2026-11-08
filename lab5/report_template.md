# Evolutionary Computation - Assignment 5 (list of improving moves in local search)

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

Neighborhood: inter-route node exchange + intra-route two-edges exchange, steepest, random start.

### Baseline: steepest local search

```
x = random solution
repeat:
    evaluate all moves of x
    if best delta < 0: apply best move
until no improving move
```

### Moves

```
node exchange (v, u, p, nx):   v in x between p and nx, u not in x
    delta = d(p,u) + d(u,nx) - d(p,v) - d(v,nx) + c(u) - c(v)

edges exchange (a -> an, b -> bn):   remove a -> an, b -> bn, add (a,b), (an,bn)
    delta = d(a,b) + d(an,bn) - d(a,an) - d(b,bn)
```

For a pair of edges the move is evaluated in two versions: (a -> an, b -> bn) and (a -> an, bn -> b).
Only one of them gives a single cycle in the current solution, the other one can become applicable
after a part of the cycle containing one of the edges is reversed.

### Steepest local search with the list of improving moves (LM)

```
x ← random solution
LM ← ∅                                        (sorted by delta, best first)
for each move m of x (edges exchange in both versions):
    if delta(m) < 0: insert m into LM

repeat
    found ← false
    for each m in LM, from the best:
        if m is node exchange (v, u, p, nx):
            if v ∈ x and u ∉ x and neighbours of v in x = {p, nx}:
                found ← true
            else:
                remove m from LM
        if m is edges exchange (a → an, b → bn):
            if edge (a, an) ∉ x or edge (b, bn) ∉ x:
                remove m from LM
            else if both edges in x have the same direction as in m, or both are reversed:
                found ← true
            else:
                skip m                        (stays in LM)
        if found:
            remove m from LM
            x ← m(x)
            evaluate new moves of x, insert those with delta < 0 into LM
            break
until not found
return x
```

New moves evaluated after applying m:

```
edges exchange, new edges e1, e2:
    e1, e2 x all other edges (both versions)
    end nodes of e1, e2 x all nodes not in x
node exchange, v replaced by u between p and nx:
    (p,u), (u,nx) x all other edges (both versions)
    p, u, nx x all nodes not in x
    all nodes in x x v
```

LM is a binary heap, skipped moves are put back after browsing.

Correctness test: for every applied move real change of objective = delta, and the final solutions are
local optima of the full neighborhood. Results equal to baseline in 183/200 (TSPA) and 173/200 (TSPB)
runs, the rest differ because of ties.

Each method run 200 times from the same random solutions as in assignments 3 and 4. C++ (g++ -O2).

## Results

Cells: average (min – max). Bold = best average and best minimum in the column, underlined = methods
from the current assignment.

### Objective function

Methods from assignments 1-4: values from assignments 3 and 4.

{{OBJ_TABLE}}

### Running time [ms] per run

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

lm.cpp (C++, standard library only). Link: https://github.com/tomasz-pres/evolutionary-computing.git

## Conclusions

- LM is about 2x faster than the baseline (TSPA 3.57 -> 1.74 ms, TSPB 3.35 -> 1.70 ms), same quality
  (differences only because of ties)
- 16x fewer evaluated moves (~124k vs ~2M), but the time gain is smaller because of the heap
  operations and checking moves that are no longer valid
- candidate moves (assignment 4) were faster (0.63 ms) but worse on TSPA, LM does not lose quality
