# Evolutionary Computation - Assignment 6 (MSLS and ILS)

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

Local search (LS) in both methods: steepest local search with the list of improving moves from
assignment 5 (inter-route node exchange + intra-route two-edges exchange).

### Multiple start local search (MSLS)

```
best ← none
repeat 400 times:
    x ← random solution
    x ← LS(x)
    if f(x) < f(best): best ← x
return best
```

### Iterated local search (ILS)

```
x ← random solution
x ← LS(x)
repeat until time limit (average time of MSLS on this instance):
    y ← perturb(x)
    y ← LS(y)
    if f(y) < f(x): x ← y
return x
```

### Perturbation

```
perturb(x):
    y ← copy of x
    # 1 double bridge
    choose random positions 0 < i < j < k < K
    cut y into 4 parts A = y[0..i-1], B = y[i..j-1], C = y[j..k-1], D = y[k..K-1]
    y ← A + C + B + D
    # 8 node replacements
    repeat 8 times:
        v ← random node in y
        u ← random node not in y
        replace v by u (u takes the position of v)
    return y
```

Double bridge changes 3 edges in a way that two-edges exchange cannot easily undo, node replacements
change the set of selected nodes. Strength was chosen in a short test (5 runs, different seed) of
0-2 double bridges and 0-16 replacements: 1 double bridge + 8 replacements was the best on average on
both instances, but the differences between 1-2 double bridges and 4-16 replacements were small.

Each method run 20 times on each instance. MSLS: 400 LS runs. ILS: stops after the average running
time of MSLS on the same instance. C++ (g++ -O2).

## Results

Cells: average (min – max). Bold = best average and best minimum in the column, underlined = methods
from the current assignment.

### Objective function

Best constructive heuristic and best local search versions from assignments 2, 3 and 5 (200 runs),
MSLS and ILS (20 runs).

{{OBJ_TABLE}}

### Running time [ms] per run

{{TIME_TABLE}}

### Number of LS runs in ILS

{{LSRUNS_TABLE}}

## Best solutions

Color of a node = its cost (bigger, outlined dots are the selected nodes).

{{GRIDS}}

## Checker

Yes, the best solutions were checked with the solution checker (in C++ and independently in Python from
the coordinates), values are the same.

## Code

ils.cpp (C++, standard library only). Link: https://github.com/tomasz-pres/evolutionary-computing.git

## Conclusions

- ILS is clearly better than MSLS in the same time: TSPA 71607 vs 73614 (-2.7%), TSPB 48255 vs 50140
  (-3.8%), the worst ILS run is better than the best MSLS run on both instances
- ILS does ~3800 LS runs in the time of 400 in MSLS, because LS after a small perturbation starts from
  an almost good solution and needs only a few moves
- MSLS is better than a single LS (best of 400 random starts), but every LS starts from scratch, so it
  can not get much better with more iterations
- ILS is the best method so far, better than all constructive heuristics and single LS versions
